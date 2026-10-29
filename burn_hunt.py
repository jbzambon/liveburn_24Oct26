#!/usr/bin/env python
"""Hunt for a satellite / radar signature of the 2026-10-24 Raleigh live-burn.

Site: 4325 Inwood Rd, Raleigh NC (35.723114 N, 78.673519 W).
Full burndown 17:00-18:00 UTC (13-14 EDT); analysis window 16:00-19:30 UTC.

Usage:  python burn_hunt.py [goes|fdc|radar|firms|viirs|tdwr|all]
Outputs go to ./burn_out/.
"""
import os
import sys
import io
import re
import json
import glob
import datetime as dt
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import s3fs
import h5netcdf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from goesgeo import latlon_to_xy, xy_to_latlon

LAT, LON = 35.723114, -78.673519
SITE_GROUND_M = 103.7   # ground elevation at the site, m MSL (USGS 3DEP)
DAY = dt.date(2026, 10, 24)
DOY = 297
T0 = pd.Timestamp("2026-10-24 16:00", tz="UTC")
T1 = pd.Timestamp("2026-10-24 19:30", tz="UTC")
BURN0 = pd.Timestamp("2026-10-24 17:00", tz="UTC")
BURN1 = pd.Timestamp("2026-10-24 18:00", tz="UTC")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "burn_out")
CACHE = os.environ.get("BURN_CACHE", os.path.join(OUT, "_cache"))
os.makedirs(OUT, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)

FS = s3fs.S3FileSystem(anon=True, config_kwargs={"read_timeout": 30, "connect_timeout": 10,
                                                  "retries": {"max_attempts": 6}})


def edt(ts):
    return (pd.Timestamp(ts) - pd.Timedelta(hours=4)).strftime("%H:%M EDT")


def shade_burn(ax):
    ax.axvspan(BURN0, BURN1, color="orange", alpha=0.15, label="burndown 17-18Z")


# ----------------------------------------------------------------------------
# GOES-19 ABI
# ----------------------------------------------------------------------------
GOES_HALF = {7: 15, 14: 15, 2: 60}   # half-width of subset in pixels per band


def goes_start(path):
    m = re.search(r"_s(\d{4})(\d{3})(\d{2})(\d{2})(\d{2})(\d)", path)
    y, j, H, M, S, _ = m.groups()
    return pd.Timestamp(f"{y}-01-01", tz="UTC") + pd.Timedelta(days=int(j) - 1, hours=int(H),
                                                              minutes=int(M), seconds=int(S))


def goes_list(product, band=None, sector=None):
    files = []
    for h in range(T0.hour, T1.hour + 1):
        d = f"noaa-goes19/{product}/2026/{DOY:03d}/{h:02d}"
        try:
            names = FS.ls(d)
        except FileNotFoundError:
            continue
        for f in names:
            b = os.path.basename(f)
            if band is not None and f"C{band:02d}_" not in b:
                continue
            if sector is not None and f"-{sector}-" not in b and f"{sector}-M6" not in b:
                continue
            t = goes_start(b)
            if T0 <= t <= T1:
                files.append((t, f))
    return sorted(files)


def _unsigned(a, attrs):
    if str(attrs.get("_Unsigned", "")).lower() == "true" and a.dtype.kind == "i":
        a = a.view(a.dtype.str.replace("i", "u"))
    return a


def _scaled(nc, name, sl=None):
    v = nc[name]
    raw = v[sl] if sl is not None else v[...]
    raw = _unsigned(np.asarray(raw), v.attrs)
    out = raw.astype("f8")
    fill = v.attrs.get("_FillValue")
    if fill is not None:
        fill = _unsigned(np.atleast_1d(np.asarray(fill)), v.attrs)[0]
        out[raw == fill] = np.nan
    out = out * float(np.ravel(v.attrs.get("scale_factor", 1.0))[0]) + \
        float(np.ravel(v.attrs.get("add_offset", 0.0))[0])
    return out


def _proj(nc):
    a = nc["goes_imager_projection"].attrs
    return {k: float(np.ravel(a[k])[0]) for k in
            ["semi_major_axis", "semi_minor_axis", "perspective_point_height",
             "longitude_of_projection_origin"]}


def goes_subset(path, half, var="Rad"):
    """Read a small window around the site over S3 (byte-range reads only)."""
    with FS.open(path, block_size=2**18, cache_type="readahead") as fo, \
            h5netcdf.File(fo, "r") as nc:
        p = _proj(nc)
        x = _scaled(nc, "x")
        y = _scaled(nc, "y")
        sx, sy = latlon_to_xy(LAT, LON, p)
        ix = int(np.argmin(np.abs(x - sx)))
        iy = int(np.argmin(np.abs(y - sy)))
        if not (half <= ix < len(x) - half and half <= iy < len(y) - half):
            return None
        sl = (slice(iy - half, iy + half + 1), slice(ix - half, ix + half + 1))
        out = {"t": goes_start(path), "file": os.path.basename(path)}
        data = _scaled(nc, var, sl)
        if var == "Rad":
            band = int(np.ravel(nc["band_id"][...])[0])
            if band <= 6:
                data = data * float(np.ravel(nc["kappa0"][...])[0])   # reflectance factor
            else:
                fk1, fk2, bc1, bc2 = (float(np.ravel(nc[k][...])[0]) for k in
                                      ["planck_fk1", "planck_fk2", "planck_bc1", "planck_bc2"])
                with np.errstate(invalid="ignore", divide="ignore"):
                    data = (fk2 / np.log(fk1 / data + 1) - bc1) / bc2
            out["dqf"] = np.asarray(nc["DQF"][sl])
        else:
            for extra in ("Power", "Temp", "Area"):
                if extra in nc.variables:
                    out[extra] = _scaled(nc, extra, sl)
        out["data"] = data
        lat, lon = xy_to_latlon(x[sl[1]], y[sl[0]], p)
        out["lat"], out["lon"] = lat, lon
        out["site_offset"] = ((sx - x[ix]) / (x[1] - x[0]), (sy - y[iy]) / (y[1] - y[0]))
        return out


def dist_km(lat, lon, lat0=LAT, lon0=LON):
    r = 6371.0
    p1, p2 = np.deg2rad(lat0), np.deg2rad(lat)
    dphi = p2 - p1
    dl = np.deg2rad(lon - lon0)
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(a))


def ring_stats(field, inner=3, outer=8):
    n = field.shape[0]
    c = n // 2
    jj, ii = np.mgrid[0:n, 0:n]
    cheb = np.maximum(np.abs(jj - c), np.abs(ii - c))
    ring = field[(cheb >= inner) & (cheb <= outer)]
    ring = ring[np.isfinite(ring)]
    med = np.median(ring)
    mad = 1.4826 * np.median(np.abs(ring - med))
    return med, mad, np.percentile(ring, 99)


def fetch_parallel(items, fn, workers=16):
    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(fn, items))


def fetch_procs(items, fn, workers=6):
    # h5netcdf/h5py reading through s3fs file objects deadlocks under threads,
    # so subset reads run in spawned processes, each with its own S3 client.
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("spawn")) as ex:
        return list(ex.map(fn, items, chunksize=2))


def _goes_triplet(args):
    tag, t, f7, f14, f2 = args
    try:
        r7 = goes_subset(f7, GOES_HALF[7])
        r14 = goes_subset(f14, GOES_HALF[14])
        r2 = goes_subset(f2, GOES_HALF[2]) if f2 else None
        if t.minute % 10 == 0:
            print(f"  [goes {tag}] {t:%H:%M}", flush=True)
        return {"t": t, "b7": r7["data"], "b14": r14["data"], "dqf7": r7["dqf"],
                "b2": None if r2 is None else r2["data"],
                "lat": r7["lat"], "lon": r7["lon"], "off": r7["site_offset"],
                "file7": r7["file"]}
    except Exception as e:  # noqa
        print("  fail", t, e, flush=True)
        return None


def goes_series(product, sector, tag):
    """Band 7/14/2 time series at the site vs a background ring."""
    cache = os.path.join(CACHE, f"goes_{tag}.npz")
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return z["recs"].tolist()
    lists = {b: dict(goes_list(product, b, sector)) for b in (7, 14, 2)}
    times = sorted(set(lists[7]) & set(lists[14]))
    print(f"[goes {tag}] {len(times)} scan times, band2 for {len(set(times) & set(lists[2]))}")

    jobs = [(tag, t, lists[7][t], lists[14][t], lists[2].get(t)) for t in times]
    recs = [r for r in fetch_procs(jobs, _goes_triplet) if r is not None]
    np.savez_compressed(cache, recs=np.array(recs, dtype=object))
    return recs


def goes_table(recs):
    rows = []
    for r in recs:
        b7, b14 = r["b7"], r["b14"]
        d = b7 - b14
        c = b7.shape[0] // 2
        m7, s7, _ = ring_stats(b7)
        m14, s14, _ = ring_stats(b14)
        md, sd, p99d = ring_stats(d)
        box = slice(c - 1, c + 2)
        row = dict(time_utc=r["t"], time_edt=edt(r["t"]),
                   bt7_site=b7[c, c], bt7_site3x3max=np.nanmax(b7[box, box]),
                   bt7_ring_med=m7, bt7_ring_mad=s7,
                   bt14_site=b14[c, c], bt14_ring_med=m14, bt14_ring_mad=s14,
                   d_site=d[c, c], d_site3x3max=np.nanmax(d[box, box]),
                   d_ring_med=md, d_ring_mad=sd, d_ring_p99=p99d,
                   dqf7_site=int(r["dqf7"][c, c]))
        row["ex7"] = row["bt7_site"] - m7
        row["ex7_3x3"] = row["bt7_site3x3max"] - m7
        row["ex14"] = row["bt14_site"] - m14
        row["exD"] = row["d_site"] - md
        row["exD_3x3"] = row["d_site3x3max"] - md
        row["zD"] = row["exD"] / max(sd, 0.05)
        if r["b2"] is not None:
            b2 = r["b2"]
            c2 = b2.shape[0] // 2
            row["refl2_site"] = np.nanmean(b2[c2 - 2:c2 + 2, c2 - 2:c2 + 2])
            row["refl2_med30km"] = np.nanmedian(b2)
            row["refl2_frac_gt0.3"] = np.mean(b2 > 0.3)
        rows.append(row)
    df = pd.DataFrame(rows).sort_values("time_utc")
    # cloud screen: band-14 colder than the clear-sky value or bright band-2
    clear14 = df["bt14_site"].quantile(0.9)
    df["cloudy_site"] = (df["bt14_site"] < clear14 - 4) | (df.get("refl2_site", 0) > 0.25)
    return df


def plot_goes_ts(df, tag, res):
    fig, ax = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
    t = df["time_utc"]
    ax[0].plot(t, df["bt7_site"], label="BT7 site px")
    ax[0].plot(t, df["bt7_site3x3max"], label="BT7 site 3x3 max", lw=0.8)
    ax[0].plot(t, df["bt7_ring_med"], label="BT7 ring median (3-8 px)")
    ax[0].fill_between(t, df["bt7_ring_med"] - 2 * df["bt7_ring_mad"],
                       df["bt7_ring_med"] + 2 * df["bt7_ring_mad"], alpha=0.2)
    ax[0].plot(t, df["bt14_site"], label="BT14 site", ls="--")
    ax[0].plot(t, df["bt14_ring_med"], label="BT14 ring median", ls="--")
    ax[0].set_ylabel("K")
    ax[0].legend(fontsize=7, ncol=3)
    ax[1].plot(t, df["ex7"], label="BT7 site - ring")
    ax[1].plot(t, df["ex14"], label="BT14 site - ring")
    ax[1].plot(t, df["exD"], label="(BT7-BT14) site - ring", lw=2)
    ax[1].plot(t, df["exD_3x3"], label="(BT7-BT14) 3x3 max - ring", lw=0.8)
    ax[1].axhline(0, color="k", lw=0.5)
    ax[1].set_ylabel("K excess")
    ax[1].legend(fontsize=7, ncol=2)
    ax[2].plot(t, df["zD"], label="z = exD / ring MAD")
    ax[2].axhline(3, color="r", ls=":")
    ax[2].set_ylabel("robust z")
    ax[2].legend(fontsize=7)
    if "refl2_site" in df:
        ax[3].plot(t, df["refl2_site"], label="0.64um refl site (2km)")
        ax[3].plot(t, df["refl2_med30km"], label="0.64um refl median +/-30km")
        ax[3].plot(t, df["refl2_frac_gt0.3"], label="frac >0.3 (cloud cover proxy)")
    ax[3].scatter(t[df["cloudy_site"]], np.full(df["cloudy_site"].sum(), 0.02), marker="|",
                  color="k", label="site flagged cloudy")
    ax[3].legend(fontsize=7)
    for a in ax:
        shade_burn(a)
        a.grid(alpha=0.3)
    ax[3].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax[0].set_title(f"GOES-19 ABI {tag} ({res}) at 4325 Inwood Rd, 2026-10-24")
    ax[3].set_xlabel("UTC")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"goes_{tag}_timeseries.png"), dpi=110)
    plt.close(fig)


def plot_goes_mosaic(recs, tag, step_min=15):
    sel = []
    last = None
    for r in sorted(recs, key=lambda r: r["t"]):
        if r["t"] < pd.Timestamp("2026-10-24 16:30", tz="UTC"):
            continue
        if last is None or (r["t"] - last) >= pd.Timedelta(minutes=step_min - 0.5):
            sel.append(r)
            last = r["t"]
    n = len(sel)
    nc = 6
    nr = int(np.ceil(n / nc))
    for kind in ("b7anom", "b2"):
        fig, axs = plt.subplots(nr, nc, figsize=(2.3 * nc, 2.4 * nr))
        for a in axs.flat:
            a.axis("off")
        for a, r in zip(axs.flat, sel):
            if kind == "b7anom":
                f = r["b7"] - np.nanmedian(r["b7"])
                im = a.imshow(f, cmap="RdBu_r", vmin=-6, vmax=6)
                c = f.shape[0] // 2
            else:
                if r["b2"] is None:
                    continue
                f = r["b2"]
                im = a.imshow(f, cmap="gray", vmin=0, vmax=0.8)
                c = f.shape[0] // 2
            a.plot(c, c, "o", mfc="none", mec="lime", ms=9, mew=1.2)
            a.set_title(f"{r['t']:%H:%M}Z", fontsize=8)
        title = ("GOES-19 band-7 anomaly (BT7 - subset median, K), ~2 km px, +/-15 px"
                 if kind == "b7anom" else
                 "GOES-19 band-2 reflectance (0.5 km px, +/-30 km)")
        fig.suptitle(f"{title} - {tag}; circle = site", fontsize=10)
        fig.colorbar(im, ax=axs, shrink=0.6)
        fig.savefig(os.path.join(OUT, f"goes_{tag}_mosaic_{kind}.png"), dpi=100)
        plt.close(fig)


def gif_goes(recs, tag):
    import imageio.v2 as imageio
    frames = []
    for r in sorted(recs, key=lambda r: r["t"]):
        if not (pd.Timestamp("2026-10-24 16:30", tz="UTC") <= r["t"]
                <= pd.Timestamp("2026-10-24 19:00", tz="UTC")):
            continue
        fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.6))
        b7 = r["b7"] - np.nanmedian(r["b7"])
        d = (r["b7"] - r["b14"])
        d = d - np.nanmedian(d)
        c = b7.shape[0] // 2
        ax[0].imshow(b7, cmap="RdBu_r", vmin=-6, vmax=6)
        ax[0].set_title("BT7 anomaly (K)")
        ax[1].imshow(d, cmap="RdBu_r", vmin=-4, vmax=4)
        ax[1].set_title("(BT7-BT14) anomaly (K)")
        if r["b2"] is not None:
            ax[2].imshow(r["b2"], cmap="gray", vmin=0, vmax=0.8)
            c2 = r["b2"].shape[0] // 2
            ax[2].plot(c2, c2, "o", mfc="none", mec="lime", ms=10)
        ax[2].set_title("0.64 um reflectance")
        for a in ax[:2]:
            a.plot(c, c, "o", mfc="none", mec="lime", ms=10)
        for a in ax:
            a.set_xticks([])
            a.set_yticks([])
        burn = BURN0 <= r["t"] < BURN1
        fig.suptitle(f"GOES-19 {tag}  {r['t']:%H:%M}Z ({edt(r['t'])})"
                     + ("  [BURNDOWN]" if burn else ""), color="red" if burn else "k")
        fig.tight_layout()
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=80)
        plt.close(fig)
        buf.seek(0)
        frames.append(imageio.imread(buf))
    imageio.mimsave(os.path.join(OUT, f"goes_{tag}_b7anom.gif"), frames, duration=0.25, loop=0)


def run_goes():
    for product, sector, tag, res, step in [("ABI-L1b-RadC", None, "conus", "5-min", 15),
                                            ("ABI-L1b-RadM", "RadM1", "meso1", "1-min", 10)]:
        recs = goes_series(product, sector, tag)
        df = goes_table(recs)
        df.to_csv(os.path.join(OUT, f"goes_{tag}_timeseries.csv"), index=False)
        plot_goes_ts(df, tag, res)
        plot_goes_mosaic(recs, tag, step)
        if tag == "meso1":
            gif_goes(recs, tag)
        print(f"[goes {tag}] wrote {len(df)} rows")


def _fdc_one(args):
    tag, t, f = args
    try:
        r = goes_subset(f, 6, var="Mask")
    except Exception as e:  # noqa
        print("  fdc fail", f, e)
        return []
    m = r["data"]
    dk = dist_km(r["lat"], r["lon"])
    c = m.shape[0] // 2
    fire = np.isin(m, [10, 11, 12, 13, 14, 15, 30, 31, 32, 33, 34, 35])
    out = [dict(product=tag, time_utc=t, kind="site_pixel", mask=int(m[c, c]),
                dist_km=float(dk[c, c]))]
    for j, i in zip(*np.where(fire & (dk <= 6))):
        out.append(dict(product=tag, time_utc=t, kind="FIRE<=6km", mask=int(m[j, i]),
                        dist_km=float(dk[j, i]), lat=float(r["lat"][j, i]),
                        lon=float(r["lon"][j, i]),
                        power_MW=float(r["Power"][j, i]), temp_K=float(r["Temp"][j, i])))
    nfar = int((fire & (dk > 6)).sum())
    if nfar:
        out.append(dict(product=tag, time_utc=t, kind="fire_beyond6km_in_13px_box", mask=nfar))
    return out


def run_fdc():
    rows = []
    for product, sector, tag in [("ABI-L2-FDCC", None, "FDCC"), ("ABI-L2-FDCM", "FDCM1", "FDCM1")]:
        files = goes_list(product, None, sector)
        if sector:
            files = [(t, f) for t, f in files if sector in f]

        for o in fetch_procs([(tag, t, f) for t, f in files], _fdc_one):
            rows.extend(o)
        print(f"[fdc] {tag}: {len(files)} files")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "goes_fdc_site.csv"), index=False)
    fires = df[df.kind == "FIRE<=6km"]
    print(f"[fdc] fire pixels within 6 km: {len(fires)}")
    print(df[df.kind == "site_pixel"].groupby(["product", "mask"]).size())
    return df


# ----------------------------------------------------------------------------
# NEXRAD KRAX Level II
# ----------------------------------------------------------------------------
def _cmap(name, fallback):
    # Py-ART >= 2.0 moved its colormaps to cmweather (names without "pyart_")
    try:
        import cmweather  # noqa: F401
    except ImportError:
        pass
    for n in (name, "pyart_" + name):
        if n in matplotlib.colormaps:
            return n
    return fallback


CM_REF = _cmap("ChaseSpectral", "turbo")
CM_RHO = _cmap("HomeyerRainbow", "viridis")
CM_ZDR = _cmap("RefDiff", "RdBu_r")
L2DIR = os.path.join(CACHE, "l2")


def l2_files():
    os.makedirs(L2DIR, exist_ok=True)
    keys = [k for k in FS.ls("unidata-nexrad-level2/2026/10/24/KRAX") if not k.endswith("_MDM")]
    out = []
    for k in keys:
        b = os.path.basename(k)
        t = pd.Timestamp(dt.datetime.strptime(b[4:19], "%Y%m%d_%H%M%S"), tz="UTC")
        if T0 - pd.Timedelta(minutes=5) <= t <= T1:
            out.append((t, k))

    def get(tk):
        dst = os.path.join(L2DIR, os.path.basename(tk[1]))
        if not os.path.exists(dst):
            FS.get(tk[1], dst)
        return tk[0], dst
    return fetch_parallel(out, get, 8)


def site_xy(radar):
    import pyart
    x, y = pyart.core.geographic_to_cartesian_aeqd(
        np.array([LON]), np.array([LAT]),
        radar.longitude["data"][0], radar.latitude["data"][0])
    return float(x[0]), float(y[0])


# control points at the same range from KRAX as the site, other azimuths
CONTROL_AZ_OFFSETS = [-60, -45, -30, 30, 45, 60, 90, 180]


MAXR = 40000.0   # only gates within 40 km of KRAX (site is at ~18 km); keeps memory small


def ngates(radar):
    return int(np.searchsorted(radar.range["data"], MAXR))


def gate_xyz(radar, s, edges=False):
    """Gate x, y (m from KRAX) and z as height AGL at the site (Py-ART z is above the antenna)."""
    n = ngates(radar) + (1 if edges else 0)
    x, y, z = [a[:, :n] for a in radar.get_gate_x_y_z(s, edges=edges)]
    return x, y, z + float(radar.altitude["data"][0]) - SITE_GROUND_M


def sweep_fields(radar, s):
    sl = radar.get_slice(s)
    n = ngates(radar)

    def f(name):
        if name not in radar.fields:
            return None
        return np.ma.filled(radar.fields[name]["data"][sl, :n].astype("f4"), np.nan)
    return f("reflectivity"), f("cross_correlation_ratio"), f("differential_reflectivity")


def gate_stats(ref, rho, zdr, mask, prefix):
    r = ref[mask]
    echo = np.isfinite(r) & (r >= -5)
    plume = echo & (r <= 25) & np.isfinite(rho[mask]) & (rho[mask] < 0.7)
    d = {f"{prefix}_ngates": int(mask.sum()), f"{prefix}_necho": int(echo.sum()),
         f"{prefix}_nplume": int(plume.sum()),
         f"{prefix}_refmax": float(np.nanmax(r)) if echo.any() else np.nan,
         f"{prefix}_refmed": float(np.nanmedian(r[echo])) if echo.any() else np.nan,
         f"{prefix}_rhomed": float(np.nanmedian(rho[mask][echo])) if echo.any() else np.nan,
         f"{prefix}_zdrstd": float(np.nanstd(zdr[mask][echo])) if echo.sum() > 2 else np.nan}
    return d


WEDGES = {"down225": 225, "up045": 45, "x135": 135, "x315": 315}


def unique_sweeps(radar):
    """First sweep at each elevation (surveillance cut, which carries dual-pol + REF)."""
    seen, out = set(), []
    for s, e in enumerate(radar.fixed_angle["data"]):
        k = round(float(e), 1)
        if k not in seen:
            seen.add(k)
            out.append((s, k))
    return out


def run_radar():
    import pyart
    import imageio.v2 as imageio
    vols = l2_files()
    print(f"[radar] {len(vols)} KRAX volumes")
    rows, tilt_rows, frames = [], [], []
    baseline = []   # lowest-tilt echo maps on a fixed polar grid, pre-burn
    grids = []
    for t, path in sorted(vols):
        radar = pyart.io.read_nexrad_archive(path)
        sx, sy = site_xy(radar)
        rng_site = np.hypot(sx, sy)
        az_site = np.rad2deg(np.arctan2(sx, sy)) % 360
        vcp = radar.metadata.get("vcp_pattern", "")
        sweeps = unique_sweeps(radar)
        for s, elev in sweeps:
            gx, gy, gz = gate_xyz(radar, s)
            ref, rho, zdr = sweep_fields(radar, s)
            if ref is None or rho is None:
                continue
            dd = np.hypot(gx - sx, gy - sy)
            row = dict(time_utc=t, time_edt=edt(t), vcp=vcp, elev=elev,
                       beam_agl_site_m=float(np.nanmedian(gz[dd < 1000])) if (dd < 1000).any() else np.nan)
            row.update(gate_stats(ref, rho, zdr, dd <= 2000, "r2km"))
            row.update(gate_stats(ref, rho, zdr, dd <= 6000, "r6km"))
            # controls: same range, other azimuths, 2 km radius
            cn = []
            for off in CONTROL_AZ_OFFSETS:
                a = np.deg2rad(az_site + off)
                cx, cy = rng_site * np.sin(a), rng_site * np.cos(a)
                cm = np.hypot(gx - cx, gy - cy) <= 2000
                cn.append(gate_stats(ref, rho, zdr, cm, "c")["c_nplume"])
            row["ctrl2km_nplume_mean"] = float(np.mean(cn))
            row["ctrl2km_nplume_max"] = float(np.max(cn))
            # wedge test: plume should appear downwind (KRDU wind from 030-060 deg -> toward SW)
            brg = np.rad2deg(np.arctan2(gx - sx, gy - sy)) % 360
            ring = (dd >= 500) & (dd <= 8000)
            for wname, wc in WEDGES.items():
                wm = ring & (np.abs((brg - wc + 180) % 360 - 180) <= 30)
                r_, h_ = ref[wm], rho[wm]
                row[f"w{wname}_n15"] = int(np.sum(r_ >= 15))
                row[f"w{wname}_n15_lowrho"] = int(np.sum((r_ >= 15) & (h_ < 0.8)))
                e_ = np.isfinite(r_) & (r_ >= 5)
                row[f"w{wname}_rhomed5"] = float(np.nanmedian(h_[e_])) if e_.any() else np.nan
            tilt_rows.append(row)
            if s == sweeps[0][0]:
                rows.append(row)
                # polar-grid echo map (0.5 deg x 250 m) for persistent-clutter screening
                az = radar.get_azimuth(s)
                ai = (np.round(az * 2).astype(int) % 720)
                g = np.full((720, ref.shape[1]), np.nan, dtype="f4")
                g[ai] = ref
                grids.append((t, g, gx.astype("f4"), gy.astype("f4"), ref, rho, zdr, ai))
        # zoomed PPIs and GIF frames of the lowest tilt
        s0 = sweeps[0][0]
        gxe, gye, _ = gate_xyz(radar, s0, edges=True)
        ref, rho, zdr = sweep_fields(radar, s0)
        X, Y = (gxe - sx) / 1000, (gye - sy) / 1000
        if pd.Timestamp("2026-10-24 16:30", tz="UTC") <= t <= pd.Timestamp("2026-10-24 19:05", tz="UTC"):
            # 4th panel: 1.3 deg REF (~450 m AGL at the site), where the plume is strongest
            s13 = min(sweeps, key=lambda se: abs(se[1] - 1.3))[0]
            gxu, gyu, _ = gate_xyz(radar, s13, edges=True)
            ref13 = sweep_fields(radar, s13)[0]
            fig, ax = plt.subplots(1, 4, figsize=(20, 5))
            for a, fld, cm, lo, hi, lab, XX, YY in [
                    (ax[0], ref, CM_REF, -10, 40, f"REF dBZ {sweeps[0][1]:.1f} deg", X, Y),
                    (ax[1], rho, CM_RHO, 0.2, 1.0, f"RHOHV {sweeps[0][1]:.1f} deg", X, Y),
                    (ax[2], zdr, CM_ZDR, -4, 8, f"ZDR dB {sweeps[0][1]:.1f} deg", X, Y),
                    (ax[3], ref13, CM_REF, -10, 40, "REF dBZ 1.3 deg",
                     (gxu - sx) / 1000, (gyu - sy) / 1000)]:
                pm = a.pcolormesh(XX, YY, np.ma.masked_invalid(fld), cmap=cm, vmin=lo, vmax=hi)
                a.plot(0, 0, "k*", ms=14, mfc="none", mew=1.5)
                for rr in (2, 6):
                    a.add_patch(plt.Circle((0, 0), rr, fill=False, ls=":", color="k"))
                a.set_xlim(-15, 15)
                a.set_ylim(-15, 15)
                a.set_aspect("equal")
                a.set_title(lab)
                a.set_xlabel("km E of site")
                fig.colorbar(pm, ax=a, shrink=0.75)
            ax[0].set_ylabel("km N of site")
            burn = BURN0 <= t < BURN1
            fig.suptitle(f"KRAX {sweeps[0][1]:.1f} deg  {t:%H:%M:%S}Z ({edt(t)})  site range "
                         f"{rng_site/1000:.1f} km az {az_site:.0f}  beam ~"
                         f"{rows[-1]['beam_agl_site_m']:.0f} m AGL" + ("  [BURNDOWN]" if burn else ""),
                         color="red" if burn else "k")
            fig.tight_layout()
            if pd.Timestamp("2026-10-24 17:00", tz="UTC") <= t <= pd.Timestamp("2026-10-24 19:00", tz="UTC"):
                fig.savefig(os.path.join(OUT, f"krax_ppi_{t:%H%M%S}.png"), dpi=90)
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=70)
            buf.seek(0)
            frames.append(imageio.imread(buf))
            plt.close(fig)
            # upper tilts: REF small multiples
            n = len(sweeps)
            fig, axs = plt.subplots(2, int(np.ceil(n / 2)), figsize=(3 * np.ceil(n / 2), 6.2))
            for a in axs.flat:
                a.axis("off")
            for a, (s, e) in zip(axs.flat, sweeps):
                gxe2, gye2, _ = gate_xyz(radar, s, edges=True)
                rf, rh, _ = sweep_fields(radar, s)
                a.axis("on")
                a.pcolormesh((gxe2 - sx) / 1000, (gye2 - sy) / 1000, np.ma.masked_invalid(rf),
                             cmap=CM_REF, vmin=-10, vmax=40)
                a.plot(0, 0, "k*", ms=10, mfc="none")
                a.set_xlim(-15, 15)
                a.set_ylim(-15, 15)
                a.set_aspect("equal")
                a.set_title(f"{e:.1f} deg", fontsize=8)
                a.tick_params(labelsize=6)
            fig.suptitle(f"KRAX REF all tilts {t:%H:%M}Z, +/-15 km around site")
            fig.tight_layout()
            if pd.Timestamp("2026-10-24 17:00", tz="UTC") <= t <= pd.Timestamp("2026-10-24 19:00", tz="UTC"):
                fig.savefig(os.path.join(OUT, f"krax_tilts_{t:%H%M%S}.png"), dpi=80)
            plt.close(fig)
        print(f"  {t:%H:%M:%S} vcp={vcp} sweeps={len(sweeps)} 2km echo={rows[-1]['r2km_necho']} "
              f"plume-like={rows[-1]['r2km_nplume']} 6km plume-like={rows[-1]['r6km_nplume']}")
    if frames:
        imageio.mimsave(os.path.join(OUT, "krax_plume.gif"), frames, duration=0.8, loop=0)

    # persistent-echo screen: gates with echo in >=50% of pre-burn (16:00-16:55Z) volumes
    pre = [g for (t, g, *_) in grids if t < pd.Timestamp("2026-10-24 16:55", tz="UTC")]
    ng = min(g.shape[1] for g in [x[1] for x in grids])
    freq = np.mean([np.isfinite(g[:, :ng]) & (g[:, :ng] >= -5) for g in pre], axis=0)
    persistent = freq >= 0.5
    sx, sy = site_xy(pyart.io.read_nexrad_archive(sorted(vols)[0][1]))
    for row, (t, g, gx, gy, ref, rho, zdr, ai) in zip(rows, grids):
        pers = persistent[ai][:, :ref.shape[1]] if ref.shape[1] <= ng else np.pad(
            persistent[ai], ((0, 0), (0, ref.shape[1] - ng)))
        new = ~pers & np.isfinite(ref) & (ref >= -5) & (ref <= 25) & np.isfinite(rho) & (rho < 0.7)
        dd = np.hypot(gx - sx, gy - sy)
        row["r2km_nplume_nonpersist"] = int((new & (dd <= 2000)).sum())
        row["r6km_nplume_nonpersist"] = int((new & (dd <= 6000)).sum())
        row["r6km_nonpersist_echo_any"] = int((~pers & np.isfinite(ref) & (ref >= -5) & (dd <= 6000)).sum())

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "krax_lowest_tilt_site.csv"), index=False)
    pd.DataFrame(tilt_rows).to_csv(os.path.join(OUT, "krax_all_tilts_site.csv"), index=False)
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    ax[0].plot(df.time_utc, df.r2km_nplume, "o-", label="plume-like gates <=2 km")
    ax[0].plot(df.time_utc, df.r6km_nplume, "o-", label="plume-like gates <=6 km")
    ax[0].plot(df.time_utc, df.r6km_nplume_nonpersist, "s--", label="<=6 km, excl. persistent echo")
    ax[0].plot(df.time_utc, df.ctrl2km_nplume_mean, "x:", label="control pts (same range) 2 km mean")
    ax[0].plot(df.time_utc, df.ctrl2km_nplume_max, "x:", label="control 2 km max")
    ax[0].set_ylabel("gates (REF -5..25, RHOHV<0.7)")
    ax[0].legend(fontsize=7)
    ax[1].plot(df.time_utc, df.r2km_refmax, "o-", label="REF max <=2 km")
    ax[1].plot(df.time_utc, df.r6km_refmed, "o-", label="REF median echo <=6 km")
    ax[1].set_ylabel("dBZ")
    ax[1].legend(fontsize=7)
    ax[2].plot(df.time_utc, df.r2km_rhomed, "o-", label="RHOHV median echo <=2 km")
    ax[2].plot(df.time_utc, df.r6km_rhomed, "o-", label="RHOHV median echo <=6 km")
    ax2 = ax[2].twinx()
    ax2.plot(df.time_utc, df.r6km_zdrstd, "k^:", label="ZDR std <=6 km")
    ax2.set_ylabel("ZDR std dB")
    ax[2].legend(fontsize=7, loc="upper left")
    for a in ax:
        shade_burn(a)
        a.grid(alpha=0.3)
    ax[2].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax[0].set_title("KRAX lowest tilt near 4325 Inwood Rd")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "krax_timeseries.png"), dpi=110)
    plt.close(fig)
    # wedge comparison, lowest tilt + all tilts for the downwind wedge
    td = pd.DataFrame(tilt_rows)
    fig, ax = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    for w in WEDGES:
        ax[0].plot(df.time_utc, df[f"w{w}_n15"], "o-", label=f"wedge {w}")
        ax[1].plot(df.time_utc, df[f"w{w}_rhomed5"], "o-", label=f"wedge {w}")
    ax[0].set_ylabel("gates REF>=15 dBZ\n(0.5-8 km, +/-30 deg)")
    ax[1].set_ylabel("median RHOHV\n(REF>=5)")
    for e, g in td.groupby("elev"):
        if e <= 4.1:
            ax[2].plot(g.time_utc, g["wdown225_n15"], ".-",
                       label=f"{e:.1f} deg (~{g.beam_agl_site_m.median():.0f} m AGL at site)")
    ax[2].set_ylabel("downwind gates REF>=15")
    for a in ax:
        shade_burn(a)
        a.grid(alpha=0.3)
        a.legend(fontsize=7, ncol=2)
    ax[2].xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax[0].set_title("KRAX: downwind (SW, 225) vs upwind (NE) and cross-wind wedges from the site")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "krax_wedges.png"), dpi=110)
    plt.close(fig)
    print(df[["time_utc", "elev", "r2km_necho", "r2km_nplume", "r6km_nplume",
              "r6km_nplume_nonpersist", "ctrl2km_nplume_mean", "r2km_refmax", "r6km_rhomed"]]
          .to_string())
    return df


def run_metar():
    url = ("https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py?station=RDU&data=tmpf&data=dwpf"
           "&data=drct&data=sknt&data=skyc1&data=skyl1&data=skyc2&data=skyl2&data=skyc3&data=skyl3"
           "&data=metar&year1=2026&month1=10&day1=24&year2=2026&month2=10&day2=25&tz=Etc/UTC"
           "&format=onlycomma&latlon=no&missing=M&trace=T&direct=no&report_type=3&report_type=4")
    import requests
    txt = requests.get(url, timeout=60).text
    df = pd.read_csv(io.StringIO(txt))
    df["valid"] = pd.to_datetime(df["valid"]).dt.tz_localize("UTC")
    df = df[(df.valid >= T0 - pd.Timedelta(hours=2)) & (df.valid <= T1 + pd.Timedelta(hours=1))]
    df.to_csv(os.path.join(OUT, "metar_KRDU.csv"), index=False)
    print(df[["valid", "drct", "sknt", "skyc1", "skyl1", "skyc2", "skyl2", "skyc3", "skyl3"]]
          .to_string())


# ----------------------------------------------------------------------------
# FIRMS + VIIRS / MODIS overpasses
# ----------------------------------------------------------------------------
FIRMS_FILES = {  # public NRT 7-day files (no MAP_KEY needed while 10/24 is within 7 days)
    "VIIRS_NOAA20": "noaa-20-viirs-c2/csv/J1_VIIRS_C2_USA_contiguous_and_Hawaii_7d.csv",
    "VIIRS_NOAA21": "noaa-21-viirs-c2/csv/J2_VIIRS_C2_USA_contiguous_and_Hawaii_7d.csv",
    "VIIRS_SNPP": "suomi-npp-viirs-c2/csv/SUOMI_VIIRS_C2_USA_contiguous_and_Hawaii_7d.csv",
    "MODIS": "modis-c6.1/csv/MODIS_C6_1_USA_contiguous_and_Hawaii_7d.csv",
}


def run_firms():
    import requests
    d = os.path.join(CACHE, "firms")
    os.makedirs(d, exist_ok=True)
    key = os.environ.get("FIRMS_MAP_KEY")
    frames = []
    for src, rel in FIRMS_FILES.items():
        fn = os.path.join(d, os.path.basename(rel))
        if not os.path.exists(fn):
            open(fn, "wb").write(requests.get(
                "https://firms.modaps.eosdis.nasa.gov/data/active_fire/" + rel, timeout=120).content)
        df = pd.read_csv(fn)
        df["source"] = src
        frames.append(df)
        if key:   # area API (also covers archive/SP data) when a key is available
            api = (f"https://firms.modaps.eosdis.nasa.gov/api/area/csv/{key}/{src}_NRT/"
                   f"{LON-0.15},{LAT-0.15},{LON+0.15},{LAT+0.15}/1/2026-10-24")
            txt = requests.get(api, timeout=120).text
            open(os.path.join(d, f"api_{src}.csv"), "w").write(txt)
    df = pd.concat(frames, ignore_index=True)
    df = df[df.acq_date == "2026-10-24"].copy()
    df["dist_km"] = dist_km(df.latitude.values, df.longitude.values)
    box = df[(abs(df.latitude - LAT) <= 0.15) & (abs(df.longitude - LON) <= 0.15)]
    box.to_csv(os.path.join(OUT, "firms_site_box_0.3deg.csv"), index=False)
    nc = df[(df.latitude.between(33.8, 36.6)) & (df.longitude.between(-84.3, -75.4))]
    passes = (nc.groupby(["source", "acq_time"]).agg(n=("latitude", "size"),
                                                     nearest_km=("dist_km", "min"))
              .reset_index())
    passes.to_csv(os.path.join(OUT, "firms_nc_passes_20261024.csv"), index=False)
    df.nsmallest(15, "dist_km").to_csv(os.path.join(OUT, "firms_nearest15.csv"), index=False)
    print(f"[firms] detections in 0.3 deg box on 10/24: {len(box)}")
    print(passes.to_string())
    print(df.nsmallest(5, "dist_km")[["source", "acq_time", "latitude", "longitude",
                                      "dist_km", "frp", "confidence"]].to_string())


VIIRS_PASSES = [  # from skyfield on Celestrak TLEs (see overpasses.txt)
    ("noaa-nesdis-n21-pds", "j02", "NOAA-21", "2026-10-24 17:36:20"),
    ("noaa-nesdis-snpp-pds", "npp", "S-NPP", "2026-10-24 18:13:10"),
    ("noaa-nesdis-n20-pds", "j01", "NOAA-20", "2026-10-24 18:32:50"),
    ("noaa-nesdis-n21-pds", "j02", "NOAA-21", "2026-10-24 19:16:30"),
]


def _h5get(fn, grp, name):
    import h5py
    with h5py.File(fn, "r") as h:
        v = h[f"All_Data/{grp}/{name}"][...]
        fac = h.get(f"All_Data/{grp}/{name}Factors")
        fac = None if fac is None else fac[...]
    v = v.astype("f8")
    if fac is not None:
        bad = v >= 65528
        v = v * fac[0] + fac[1]
        v[bad] = np.nan
    else:
        v[v < -999] = np.nan
    return v


def viirs_granule(bucket, sat, t):
    """Find the I-band granule(s) around time t containing the site."""
    t = pd.Timestamp(t, tz="UTC")
    vd = os.path.join(CACHE, "viirs")
    os.makedirs(vd, exist_ok=True)
    geo = FS.ls(f"{bucket}/VIIRS-IMG-GEO-TC/2026/10/24")

    def tstart(k):
        m = re.search(r"_t(\d{6})", k)
        return pd.Timestamp(f"2026-10-24 {m.group(1)[:2]}:{m.group(1)[2:4]}:{m.group(1)[4:6]}", tz="UTC")
    geo = sorted([(tstart(k), k) for k in geo if f"_{sat}_" in k])
    i0 = max(i for i, (ts, _) in enumerate(geo) if ts <= t)
    for i in (i0, i0 - 1, i0 + 1):
        ts, gk = geo[i]
        gfn = os.path.join(vd, os.path.basename(gk))
        if not os.path.exists(gfn):
            FS.get(gk, gfn)
        lat = _h5get(gfn, "VIIRS-IMG-GEO-TC_All", "Latitude")
        lon = _h5get(gfn, "VIIRS-IMG-GEO-TC_All", "Longitude")
        dk = dist_km(lat, lon)
        j, k = np.unravel_index(np.nanargmin(dk), dk.shape)
        if dk[j, k] < 1.0 and 20 <= j < lat.shape[0] - 20:
            tag = re.search(r"_d\d{8}_t\d{7}", gk).group(0)
            return gfn, tag, (j, k), lat, lon
    return None


def run_viirs():
    rows = []
    for bucket, sat, name, tov in VIIRS_PASSES:
        g = viirs_granule(bucket, sat, tov)
        if g is None:
            print(f"[viirs] {name} {tov}: site not found in granules")
            continue
        gfn, tag, (j, k), lat, lon = g
        H = 25
        win = (slice(j - H, j + H + 1), slice(k - H, k + H + 1))
        bands = {}
        for b in ("I01", "I04", "I05"):
            prod = f"VIIRS-{b.replace('I0', 'I')}-SDR"
            keys = [x for x in FS.ls(f"{bucket}/{prod}/2026/10/24") if tag in x]
            fn = os.path.join(CACHE, "viirs", os.path.basename(keys[0]))
            if not os.path.exists(fn):
                FS.get(keys[0], fn)
            grp = f"{prod}_All"
            name_ = "Reflectance" if b == "I01" else "BrightnessTemperature"
            bands[b] = _h5get(fn, grp, name_)[win]
        sza = _h5get(gfn, "VIIRS-IMG-GEO-TC_All", "SatelliteZenithAngle")[j, k]
        la, lo = lat[win], lon[win]
        dk = dist_km(la, lo)
        i4, i5, i1 = bands["I04"], bands["I05"], bands["I01"]
        d45 = i4 - i5
        near = dk <= 0.75          # site +/- ~2 pixels (geolocation & PSF)
        bg = (dk > 1.5) & (dk <= 6)
        row = dict(sat=name, granule=os.path.basename(gfn), overpass_utc=tov,
                   overpass_edt=edt(tov), sat_zenith=float(sza),
                   I4_site_max=np.nanmax(i4[near]), I4_bg_med=np.nanmedian(i4[bg]),
                   I4_bg_std=np.nanstd(i4[bg]), I5_site=np.nanmax(i5[near]),
                   I5_bg_med=np.nanmedian(i5[bg]),
                   d45_site_max=np.nanmax(d45[near]), d45_bg_med=np.nanmedian(d45[bg]),
                   d45_bg_std=np.nanstd(d45[bg]), d45_bg_p99=np.nanpercentile(d45[bg], 99),
                   I1_site=np.nanmean(i1[near]), I1_bg_med=np.nanmedian(i1[bg]))
        row["I1_bg_frac_gt0.3"] = float(np.mean(i1[bg] > 0.3))
        row["I4_excess"] = row["I4_site_max"] - row["I4_bg_med"]
        row["I4_z"] = row["I4_excess"] / row["I4_bg_std"]
        row["d45_excess"] = row["d45_site_max"] - row["d45_bg_med"]
        row["d45_z"] = row["d45_excess"] / row["d45_bg_std"]
        # rank of the site's warmest pixel among the whole +/-25 px window
        row["I4_site_rank_in_window"] = int((i4 > row["I4_site_max"]).sum()) + 1
        rows.append(row)
        fig, ax = plt.subplots(1, 4, figsize=(17, 4.2))
        for a, f, cm, lab, lo_, hi_ in [
                (ax[0], i1, "gray", "I1 0.64um refl", 0, 0.8),
                (ax[1], i4 - np.nanmedian(i4), "magma", "I4 3.74um BT anomaly (K)", -10, 30),
                (ax[2], i5 - np.nanmedian(i5), "RdBu_r", "I5 11.45um BT anomaly (K)", -8, 8),
                (ax[3], d45 - np.nanmedian(d45), "magma", "(I4-I5) anomaly (K)", -10, 30)]:
            im = a.imshow(f, cmap=cm, vmin=lo_, vmax=hi_)
            a.plot(k - win[1].start, j - win[0].start, "o", mfc="none", mec="lime", ms=12, mew=1.5)
            a.set_title(lab, fontsize=9)
            a.set_xticks([])
            a.set_yticks([])
            fig.colorbar(im, ax=a, shrink=0.8)
        fig.suptitle(f"{name} VIIRS I-bands {tov[11:16]}Z ({edt(tov)}), VZA {sza:.0f} deg, "
                     f"+/-25 px window (scan geometry, not north-up); circle = site")
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"viirs_{name}_{tov[11:13]}{tov[14:16]}.png"), dpi=100)
        plt.close(fig)
        print(f"[viirs] {name} {tov}Z VZA {sza:.0f}: I4 site max {row['I4_site_max']:.1f} K vs bg "
              f"{row['I4_bg_med']:.1f}+/-{row['I4_bg_std']:.1f}; d45 excess {row['d45_excess']:.1f} K "
              f"(z={row['d45_z']:.1f}); I1 site {row['I1_site']:.2f}")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "viirs_overpass_site.csv"), index=False)


# ----------------------------------------------------------------------------
# TDWR TRDU Level III (Unidata L3 bucket mirrors the NCEI products; Py-ART can't read
# TDWR code 180, so use MetPy)
# ----------------------------------------------------------------------------
def run_tdwr():
    import requests
    import xml.etree.ElementTree as ET
    from metpy.io import Level3File
    d = os.path.join(CACHE, "tdwr")
    os.makedirs(d, exist_ok=True)
    keys = []
    for prod in ("TZ0", "TZ1", "TZ2"):
        for h in range(16, 20):
            r = requests.get("https://unidata-nexrad-level3.s3.amazonaws.com/", timeout=60,
                             params={"list-type": "2", "prefix": f"RDU_{prod}_2026_10_24_{h:02d}"})
            ns = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}
            keys += [k.text for k in ET.fromstring(r.content).findall(".//s:Key", ns)]

    def get(k):
        fn = os.path.join(d, k)
        if not os.path.exists(fn):
            FS.get(f"unidata-nexrad-level3/{k}", fn)
        return fn
    files = fetch_parallel(keys, get, 8)
    rows, snaps = [], {}
    for fn in sorted(files):
        b = os.path.basename(fn)
        t = pd.Timestamp(dt.datetime.strptime(b[8:], "%Y_%m_%d_%H_%M_%S"), tz="UTC")
        if not (T0 <= t <= T1):
            continue
        L = Level3File(fn)
        blk = L.sym_block[0][0]
        ref = np.asarray(L.map_data(blk["data"]), dtype="f4")
        az = np.asarray(blk["start_az"]) + 0.5
        rng = (np.arange(ref.shape[1]) + 0.5) * 0.15   # km, 150 m gates
        # flat-earth offsets from TRDU are fine at ~30 km for this purpose
        A, R = np.meshgrid(np.deg2rad(az), rng, indexing="ij")
        x = R * np.sin(A)
        y = R * np.cos(A)
        sx = (LON - L.lon) * 111.32 * np.cos(np.deg2rad(LAT))
        sy = (LAT - L.lat) * 110.57
        dd = np.hypot(x - sx, y - sy)
        brg = np.rad2deg(np.arctan2(x - sx, y - sy)) % 360
        ring = (dd >= 0.5) & (dd <= 8)
        row = dict(time_utc=t, time_edt=edt(t), product=b[4:7],
                   elev=round(float(L.metadata["el_angle"]), 1))
        for wname, wc in WEDGES.items():
            wm = ring & (np.abs((brg - wc + 180) % 360 - 180) <= 30)
            row[f"w{wname}_n15"] = int(np.sum(ref[wm] >= 15))
            row[f"w{wname}_max"] = float(np.nanmax(ref[wm]))
        rows.append(row)
        if b[4:7] == "TZ0":
            snaps[t] = (x - sx, y - sy, ref)
    df = pd.DataFrame(rows).sort_values(["product", "time_utc"])
    df.to_csv(os.path.join(OUT, "tdwr_trdu_wedges.csv"), index=False)
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for (p, e), g in df.groupby(["product", "elev"]):
        ax.plot(g.time_utc, g.wdown225_n15, ".-", label=f"{p} {e} deg downwind (SW)")
        ax.plot(g.time_utc, g.wup045_n15, ":", label=f"{p} upwind (NE)")
    shade_burn(ax)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    ax.set_ylabel("gates REF>=15 dBZ (0.5-8 km, +/-30 deg)")
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
    ax.set_title("TDWR TRDU Level III: downwind vs upwind wedge from the site")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "tdwr_trdu_wedges.png"), dpi=110)
    plt.close(fig)
    sel = [t for t in sorted(snaps) if pd.Timestamp("2026-10-24 17:00", tz="UTC") <= t
           <= pd.Timestamp("2026-10-24 18:30", tz="UTC")]
    nc = 5
    nr = int(np.ceil(len(sel) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.2 * nc, 3.2 * nr))
    for a in axs.flat:
        a.axis("off")
    for a, t in zip(axs.flat, sel):
        X, Y, Z = snaps[t]
        m = (np.abs(X) < 16) & (np.abs(Y) < 16)
        a.axis("on")
        a.scatter(X[m], Y[m], c=Z[m], s=2, cmap=CM_REF, vmin=-10, vmax=40, marker="s")
        a.plot(0, 0, "k*", ms=10, mfc="none")
        a.set_xlim(-15, 15)
        a.set_ylim(-15, 15)
        a.set_aspect("equal")
        a.set_title(f"TZ0 0.3 deg {t:%H:%M}Z", fontsize=8)
        a.tick_params(labelsize=6)
    fig.suptitle("TDWR TRDU 0.3 deg reflectivity, +/-15 km around site (km)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "tdwr_trdu_tz0_montage.png"), dpi=90)
    plt.close(fig)
    print(df[df["product"] == "TZ0"][["time_utc", "wdown225_n15", "wdown225_max", "wup045_n15",
                                       "wx135_n15", "wx315_n15"]].to_string())
    print(df.pivot_table(index="time_utc", columns="product", values="wdown225_n15").to_string())


# ----------------------------------------------------------------------------
if __name__ == "__main__":
    what = sys.argv[1:] or ["all"]
    steps = {"goes": run_goes, "fdc": run_fdc, "radar": run_radar, "metar": run_metar,
             "firms": run_firms, "viirs": run_viirs, "tdwr": run_tdwr}
    for name in (list(steps) if "all" in what else what):
        steps[name]()

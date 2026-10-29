"""Dual-pol stats of KRAX plume gates in the downwind wedge (needs `burn_hunt.py radar` cache)."""
import os, sys, glob
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pyart
import burn_hunt as B

for f in sorted(glob.glob(os.path.join(B.L2DIR, "KRAX20261024_1[678]*"))):
    t = f.split("_")[-2]
    if not ("1723" <= t[:4] <= "1808" or t[:4] == "1640"):
        continue
    r = pyart.io.read_nexrad_archive(f)
    sx, sy = B.site_xy(r)
    P = {"ref": [], "rho": [], "zdr": [], "z": []}
    Q = {"rho": [], "zdr": []}
    for s, e in B.unique_sweeps(r):
        if e > 2.5:
            continue
        gx, gy, gz = B.gate_xyz(r, s)
        ref, rho, zdr = B.sweep_fields(r, s)
        dd = np.hypot(gx - sx, gy - sy) / 1000
        brg = np.rad2deg(np.arctan2(gx - sx, gy - sy)) % 360
        ring = (dd >= 0.5) & (dd <= 8)
        dw = ring & (np.abs((brg - 225 + 180) % 360 - 180) <= 30)
        m = dw & (ref >= 15)
        for k, v in zip(("ref", "rho", "zdr", "z"), (ref, rho, zdr, gz)):
            P[k] += list(v[m])
        o = ring & ~dw & (ref >= 15)
        Q["rho"] += list(rho[o])
        Q["zdr"] += list(zdr[o])
    a = {k: np.array(v) for k, v in P.items()}
    b = {k: np.array(v) for k, v in Q.items()}
    if len(a["ref"]) < 3:
        print(t, "n=", len(a["ref"]))
        continue
    print(f"{t[:2]}:{t[2:4]}Z n={len(a['ref'])} REF med {np.nanmedian(a['ref']):.1f} "
          f"RHOHV med {np.nanmedian(a['rho']):.2f} frac<0.7 {np.nanmean(a['rho'] < 0.7):.2f} "
          f"ZDR std {np.nanstd(a['zdr']):.1f} | other wedges RHOHV med {np.nanmedian(b['rho']):.2f} "
          f"ZDR std {np.nanstd(b['zdr']):.1f}")

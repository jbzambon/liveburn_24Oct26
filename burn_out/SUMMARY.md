# Live-burn signal hunt: 4325 Inwood Rd, Raleigh, 2026-10-24

Site 35.7231 N, 78.6735 W. Burndown expected 17:00–18:00 UTC (13:00–14:00 EDT).
KRDU METARs: wind from 030–060° at 10–13 kt, so the plume should drift toward the SW (~225°).
Sky: FEW/SCT before 16:50Z, then BKN at 3,800–5,000 ft through 18:51Z.

## Credible signals

| Dataset | Time UTC (EDT) | Magnitude | Notes |
|---|---|---|---|
| **KRAX Level II** (VCP 35) | first aloft 17:32Z (13:32); peak 17:40–17:58Z (13:40–13:58); back to baseline ~18:20Z (14:20) | Downwind wedge (SW, 0.5–8 km): gates ≥15 dBZ go from a baseline of 2–18 to **159** (0.2°) and **204** (1.3°). Plume gates: REF median 18–22 dBZ, max 35–40 dBZ; **RHOHV median 0.63–0.74** (background 0.87–0.93); 50–60% of gates have RHOHV < 0.7. | Elongated plume 2–9 km SW of the site, along the wind. It is absent in the upwind and cross-wind wedges. It shows at 0.2–2.4° and not at 3.1° or above, so the top is about 0.8–1 km. It first appeared at mid-levels (0.9–1.8°) at 17:32Z. ZDR is *not* distinctive: the plume's ZDR spread (std ~3 dB) is smaller than the insect background's (std ~4.5 dB). |
| **TDWR TRDU Level III** (TZ0/TZ1, 0.3°/1.0°) | 17:27–17:51Z (13:27–13:51), peak 17:45Z | Downwind gates ≥15 dBZ: baseline 0–1, up to **24**; max 25 dBZ. Upwind 0–2. | An independent radar with different viewing geometry confirms the KRAX plume. That rules out a KRAX-specific artifact. |
| **FIRMS VIIRS S-NPP** | 18:13Z (14:13) | Detection 0.18 km from the site; FRP **3.1 MW**; BT I4 330.4 K; nominal confidence | Nothing within 3 km on any other day in the 7-day file, so this is not a static hot spot. |
| **FIRMS VIIRS NOAA-20** | 18:32Z (14:32) | Detection 0.34 km from the site; FRP **1.4 MW**; BT I4 328.4 K; nominal | The decaying fire or debris pile. |
| VIIRS I-band SDR (own check) | 18:13Z and 18:32Z | The site pixel is the hottest of 2,601 in a ±25 px window. I4 is +29 K and +28 K above background (z ≈ 5–6), and I4−I5 is +24 K and +28 K. | Same pixels as FIRMS. |

## Suggestive only

* **GOES-19 M1 (1-min), 17:47Z (13:47 EDT):** at the site pixel, BT7 jumps +2.3 K for one minute while BT14 changes by +0.1 K and band-2 reflectance is unchanged (0.23). Smaller jumps occur at 17:43–17:44Z. The site pixel is the local BT7 maximum, with PSF-like falloff into its neighbours.
  * I tested every pixel-minute in the window (83,168, 16:00–19:30Z). This is the **only** one with a ≥2 K one-minute BT7 rise and flat BT14.
  * The timing matches the radar plume peak.
  * But it is a single minute under broken cloud, and I found it by looking. Treat it as possible, not confirmed. See `goes_meso1_flicker_test.txt`.

## No credible signal

* **GOES-19 CONUS (5-min) and M1 BT7 / BT7−BT14 excess:**
  * There is a positive BT7−BT14 excess during 17:00–18:26Z (+3 to +8 K, z < 2), but it is a cloud artifact. BT14 at the site drops 3–6 K and 0.64 µm reflectance rises to 0.25–0.42 over the same period, which is the usual daytime reflected-solar contamination of band 7.
  * This fails the "band-7 rise without an 11 µm change" test.
  * In the 18:30–18:55Z clear gap, the site-adjacent max BT7 reaches +2.9 K. It never exceeds the ring's 99th percentile, and it coincides with a cloud edge (BT14 falling).
* **GOES ABI-L2-FDCC / FDCM1:** no fire pixels within 6 km in any of the 42 + 210 files. The FDC algorithm itself labels the site "cloud_pixel_detected_by_albedo_threshold_test" (mask 215) for most of 17:01–18:29Z.
* **VIIRS NOAA-21, 17:36Z (13:36 EDT, view zenith 51°):** this pass fell mid-burn, while the radar plume was strong. Nothing at the site: I4 +2.4 K, z 0.6, rank 635 of 2,601. I5 at the site is warmer than the surroundings, so it was not under opaque cloud.
  * Likely reasons: at 51° the view is oblique, the 0.8 km pixel is larger than at nadir, and the flames may still have been mostly inside the structure. This is a speculative interpretation, not something the data show.
* **VIIRS NOAA-21, 19:16Z (view zenith 67°):** nothing.
* **MODIS:** no pass in the window. Terra crossed at 14:21Z and Aqua at 20:43Z. No MODIS fire within 90 km that day.
* **KRAX "within 2 km of the site" counts:** gates with REF −5 to 25 dBZ and RHOHV < 0.7 are dominated by clear-air, insect and clutter echo. Same-range control points show similar counts, so this metric alone is not useful. The downwind-wedge test is the discriminating one.

## Effect of cloud cover on the satellite results

* A cumulus deck (BKN 3.8–5 kft) built over the site from about 16:42Z and covered the whole burndown hour. It cleared only briefly around 18:30–18:55Z.
* Under that cloud, GOES band 7 is dominated by reflected sunlight from cloud tops. The fire's expected contribution is far smaller:
  * Using FRP ≈ 1.4–3 MW from VIIRS, the BT7 excess in a ~6.5 km² ABI pixel should be only about 0.3–1 K, before PSF spreading.
  * The clear-sky ring noise is ~1 K MAD.
  * So GOES could not have detected the late-phase fire even in clear air. A peak-phase fire of ~10 MW might have reached 2–3 K, but it was hidden by cloud.
* VIIRS succeeded because its 375 m pixels dilute the fire ~50× less than ABI's, and the 18:13Z and 18:32Z passes fell in a cloud gap.

## Timeline

* 17:27–17:32Z: smoke plume first seen by both radars, aloft at ~0.3–0.8 km.
* 17:40–17:58Z: plume peak.
* 17:47Z: possible single-minute GOES flicker.
* ~18:20Z: plume gone.
* 18:13Z and 18:32Z: VIIRS sees the residual hot debris, FRP 3.1 MW falling to 1.4 MW.

## Files

* `burn_hunt.py` (repo root): run as `python burn_hunt.py [goes fdc radar metar firms viirs tdwr]`.
* **GOES:** `goes_{conus,meso1}_timeseries.{csv,png}`, `goes_*_mosaic_{b7anom,b2}.png`, `goes_meso1_b7anom.gif`, `goes_fdc_site.csv`, `goes_meso1_flicker_test.txt`
* **KRAX:** `krax_lowest_tilt_site.csv`, `krax_all_tilts_site.csv`, `krax_timeseries.png`, `krax_wedges.png`, `krax_ppi_HHMMSS.png` (REF/RHOHV/ZDR at 0.2° plus REF at 1.3°), `krax_tilts_HHMMSS.png`, `krax_plume.gif`, `krax_plume_dualpol.txt`
* **TDWR:** `tdwr_trdu_wedges.{csv,png}`, `tdwr_trdu_tz0_montage.png`
* **FIRMS / VIIRS:** `firms_site_box_0.3deg.csv`, `firms_nearest15.csv`, `firms_nc_passes_20261024.csv`, `overpasses.txt`, `viirs_overpass_site.csv`, `viirs_<sat>_<HHMM>.png`
* **Context:** `metar_KRDU.csv`

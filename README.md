# liveburn_24Oct26

Satellite and radar search for the signature of a live-burn training fire (acquired-structure
burndown) at 4325 Inwood Rd, Raleigh NC (35.7231 N, 78.6735 W) on 2026-10-24, with peak
burndown 17:00-18:00 UTC.

**Findings:** see [`burn_out/SUMMARY.md`](burn_out/SUMMARY.md). In short:
- KRAX and TDWR TRDU both show a low-RHOHV smoke plume downwind (SW), 17:27-18:20Z, peaking 17:40-17:58Z.
- VIIRS (FIRMS) detected the fire at 18:13Z (S-NPP, 3.1 MW) and 18:32Z (NOAA-20, 1.4 MW).
- GOES-19 shows no credible signal because cloud covered the site (one suggestive 1-min band-7 flicker at 17:47Z).

## Running

```bash
uv venv -p 3.12 .venv && VIRTUAL_ENV=.venv uv pip install -r requirements.txt
.venv/bin/python burn_hunt.py metar goes fdc radar firms viirs tdwr
.venv/bin/python extras/flicker_test.py      # after `goes`
.venv/bin/python extras/plume_dualpol.py     # after `radar`
```

All data come from public anonymous sources: the `noaa-goes19`, `unidata-nexrad-level2`,
`unidata-nexrad-level3` and `noaa-nesdis-{n20,n21,snpp}-pds` S3 buckets, FIRMS, IEM ASOS and
Celestrak. Raw downloads are cached in `burn_out/_cache/` (about 1 GB, not committed).

Notes:
- The FIRMS public 7-day CSVs only cover a date for about a week. Set `FIRMS_MAP_KEY` to use the
  area API instead.
- GOES subset reads run in processes, not threads, because h5py over s3fs deadlocks when threaded.

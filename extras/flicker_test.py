"""GOES-19 M1 band-7 one-minute 'flicker' test at the site (needs `burn_hunt.py goes` cache)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from scipy.ndimage import median_filter
import burn_hunt as B

recs = sorted(np.load(os.path.join(B.CACHE, "goes_meso1.npz"), allow_pickle=True)["recs"].tolist(),
              key=lambda r: r["t"])
T = np.array([r["t"] for r in recs])
b7 = np.stack([r["b7"] for r in recs])
b14 = np.stack([r["b14"] for r in recs])


def local_anom(a):
    # pixel minus median of its 7x7 neighbourhood excluding the 3x3 core
    fp = np.ones((7, 7), bool)
    fp[2:5, 2:5] = False
    return a - np.stack([median_filter(x, footprint=fp, mode="nearest") for x in a])


a7, a14 = local_anom(b7), local_anom(b14)
d7, d14 = np.diff(a7, axis=0), np.diff(a14, axis=0)
c = 15
inner = (slice(None), slice(4, -4), slice(4, -4))
pop_d, pop_d14 = d7[inner].ravel(), d14[inner].ravel()
ok = np.abs(pop_d14) < 0.5
big = (pop_d >= 2) & ok
print(f"pixel-minutes with dBT7>=+2 K and |dBT14|<0.5 K: {big.sum()} of {ok.sum()}")
jj = np.where((d7[:, c, c] >= 1.5) & (np.abs(d14[:, c, c]) < 0.5))[0]
print("site minutes with dBT7>=1.5 K and |dBT14|<0.5 K:", [T[j + 1].strftime("%H:%M") for j in jj])
cnt = ((d7 >= 1.5) & (np.abs(d14) < 0.5))[:, 4:-4, 4:-4].sum(0)
print(f"site count {cnt[c-4, c-4]}; window max {cnt.max()}; rank {(cnt > cnt[c-4, c-4]).sum()+1}/{cnt.size}")

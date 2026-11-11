#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd
from scipy.stats import binomtest

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WK, "scripts"))
import panel as P

GRID = os.path.join(WK, "benchmark_grid")
PAIR = {(64, 100): ("45_spotsweeper_8um_k768", "19_spatialartifacts_8um"),
        (64, 12): ("26_spotsweeper_32um_k48", "20_spatialartifacts_16um"),
        (256, 100): ("27_spotsweeper_64um_k48", "19_spatialartifacts_8um"),
        (256, 12): ("27_spotsweeper_64um_k48", "20_spatialartifacts_16um"),
        (512, 100): ("28_spotsweeper_128um_k48", "19_spatialartifacts_8um"),
        (512, 12): ("28_spotsweeper_128um_k48", "20_spatialartifacts_16um")}

def removed(run, arm, cond):
    f = os.path.join(GRID, run, "results", arm, f"{cond}.txt")
    if not os.path.exists(f):
        return None
    return {ln.strip() for ln in open(f) if ln.strip()}

def truth(run, cond):
    d = os.path.join(GRID, run, cond)
    gt = pd.read_parquet(os.path.join(d, "ground_truth_008um.parquet"), columns=["barcode", "is_artifact"])
    uni = set(gt.barcode)
    pos = os.path.join(d, "square_008um", "spatial", "tissue_positions.parquet")
    if os.path.exists(pos):
        p = pd.read_parquet(pos)
        if "in_tissue" in p.columns:
            uni = set(p.barcode[p.in_tissue == 1])
    return set(gt.barcode[gt.is_artifact]) & uni, uni

rows = []
for sec, dep, cell in P.cells(GRID, depths=(100, 12)):
    run = os.path.basename(cell)
    for dm in (64, 256, 512):
        ss, sa = PAIR[(dm, dep)]
        conds = ["clean"] + [f"{k}{e}_{dm}um" for k in ("dropout01", "mito") for e in ("", "edge")]
        for cond in conds:
            if not os.path.isdir(os.path.join(GRID, run, cond)):
                continue
            a, b = removed(run, ss, cond), removed(run, sa, cond)
            if a is None or b is None:
                continue
            t, uni = truth(run, cond)
            a, b = a & uni, b & uni
            u = a | b
            r = dict(section=sec, chem=P.chemistry(sec), depth=dep, dm=dm, cond=cond, ss=ss, sa=sa)
            if cond == "clean":
                r.update(ss_v=len(a) / len(uni), sa_v=len(b) / len(uni), un_v=len(u) / len(uni), metric="fpr")
            else:
                r.update(ss_v=len(a & t) / len(t), sa_v=len(b & t) / len(t), un_v=len(u & t) / len(t), metric="recall")
            rows.append(r)

D = pd.DataFrame(rows)
D.to_csv(os.path.join(GRID, "figures", "union_region.csv"), index=False)

def sign(x, y):
    d = (x - y).dropna()
    up, dn = int((d > 0).sum()), int((d < 0).sum())
    return up, dn, (binomtest(up, up + dn).pvalue if up + dn else np.nan)

for ch in ("probe", "3prime"):
    for dep in (100, 12):
        for dm in (64, 256, 512):
            x = D[(D.chem == ch) & (D.depth == dep) & (D.dm == dm)]
            for cond in sorted(x.cond.unique()):
                y = x[x.cond == cond]
                s = (f"{ch:6} d{dep:<3} {dm:3} {cond:20} n={len(y):2} SS {100*y.ss_v.median():5.1f} "
                     f"SA {100*y.sa_v.median():5.1f} union {100*y.un_v.median():5.1f}")
                up, dn, p = sign(y.un_v, y.ss_v); s += f" | union>SS {up}/{up+dn} p={p:.2g}"
                up, dn, p = sign(y.un_v, y.sa_v); s += f" | union>SA {up}/{up+dn} p={p:.2g}"
                print(s)

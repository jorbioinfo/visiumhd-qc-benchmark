import argparse, glob, os, re
import numpy as np, pandas as pd
from scipy import ndimage

import panel

BC = re.compile(r"s_008um_(\d+)_(\d+)-1")
BOOT = 2000
RNG = np.random.default_rng(0)

def coords(barcodes):
    rc = np.array([[int(m.group(1)), int(m.group(2))]
                   for m in (BC.match(b) for b in barcodes) if m], dtype=np.int32)
    return rc

def regions(truth):
    rc = coords(truth)
    if not len(rc):
        return {}, 0
    r0, c0 = rc.min(0)
    grid = np.zeros((rc[:, 0].max() - r0 + 3, rc[:, 1].max() - c0 + 3), bool)
    grid[rc[:, 0] - r0 + 1, rc[:, 1] - c0 + 1] = True
    lab, n = ndimage.label(grid, structure=np.ones((3, 3)))
    ids = lab[rc[:, 0] - r0 + 1, rc[:, 1] - c0 + 1]
    return dict(zip(truth, ids)), n

def boot_ci(per_region):
    a = np.asarray(per_region, float)
    if len(a) < 2:
        return float(a.mean()) if len(a) else np.nan, np.nan, np.nan
    idx = RNG.integers(0, len(a), size=(BOOT, len(a)))
    means = a[idx].mean(1)
    return float(a.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))

def wilson(k, n, z=1.96):
    if not n:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--depths", default="100,12")
    ap.add_argument("--conds", default="dropout01_256um,mito_256um")
    ap.add_argument("--out", default="/work/benchmark_grid/figures/region_ci.csv")
    a = ap.parse_args()
    depths = [int(x) for x in a.depths.split(",")]
    conds = a.conds.split(",")
    rows = []
    for sec, dep, cell in panel.cells(a.grid, depths):
        tag = os.path.basename(cell)
        sc_path = os.path.join(cell, "scored", "scores.csv")
        sc = panel.scores(cell)
        for cond in conds:
            gt = os.path.join(cell, cond, "ground_truth_008um.parquet")
            if not os.path.exists(gt):
                continue
            g = pd.read_parquet(gt)
            truth = list(g.barcode[g.is_artifact])
            rmap, nreg = regions(truth)
            if nreg < 2:
                continue
            for arm_dir in sorted(glob.glob(os.path.join(cell, "results", "*"))):
                arm = os.path.basename(arm_dir)
                pf = os.path.join(arm_dir, f"{cond}.txt")
                if not os.path.exists(pf):
                    continue
                pred = set(open(pf).read().split())
                hit = np.zeros(nreg + 1); tot = np.zeros(nreg + 1)
                for b, rid in rmap.items():
                    tot[rid] += 1
                    if b in pred:
                        hit[rid] += 1
                keep = tot > 0
                per = hit[keep] / tot[keep]
                mean, lo, hi = boot_ci(per)
                row = sc[(sc.arm == arm) & (sc.condition == cond)]
                if not len(row):
                    continue
                rec = float(row.recall.iloc[0])
                n_total = int(row.n_total.iloc[0])
                clean = sc[(sc.arm == arm) & sc["kind"].isin(["clean", "none"])]
                if not len(clean):
                    continue
                cost = float(clean.discard_frac.iloc[0])
                clo, chi = wilson(cost * n_total, n_total)
                half_lo, half_hi = mean - lo, hi - mean
                rows.append(dict(
                    section=panel.PREFIX + sec, depth=dep, chem=panel.chemistry(sec),
                    cond=cond, arm=arm,
                    n_regions=int(keep.sum()), recall=rec,
                    recall_lo=max(0.0, rec - half_lo), recall_hi=min(1.0, rec + half_hi),
                    cost=cost, cost_lo=clo, cost_hi=chi))
        print("done", tag, flush=True)
    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    out.to_csv(a.out, index=False)
    print("wrote", a.out, len(out), "rows")

if __name__ == "__main__":
    main()

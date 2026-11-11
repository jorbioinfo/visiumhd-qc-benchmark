#!/usr/bin/env python3
import os, sys
import numpy as np, pandas as pd

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WK, "scripts"))
import panel as P

GRID = os.path.join(WK, "benchmark_grid")
DEPTHS = (100, 12)
MAD = {100: "02_mad_8um", 12: "09_mad_32um"}
SS_AP = {32: "25_spotsweeper_16um_k48", 64: "26_spotsweeper_32um_k48", 256: "27_spotsweeper_64um_k48"}
SS_WIDE = {32: "44_spotsweeper_8um_k192", 64: "45_spotsweeper_8um_k768", 256: "46_spotsweeper_8um_k3072"}
COLS = [("Section-wide (MAD)", lambda dm, dep: MAD[dep]),
        ("Neighbourhood, A&P", lambda dm, dep: SS_AP[dm]),
        ("Neighbourhood, wider", lambda dm, dep: SS_WIDE[dm]),
        ("Region, 8 µm", lambda dm, dep: "19_spatialartifacts_8um"),
        ("Region, A&P 16 µm", lambda dm, dep: "20_spatialartifacts_16um")]
ROWS = [(f"{lab}, {dm} µm{where}", k + e, dm)
        for e, where in (("", ""), ("edge", ", at the edge"))
        for k, lab in (("dropout01", "Signal loss"), ("mito", "Mito. elevation"))
        for dm in (32, 64, 256)]

def load():
    rows = []
    for sec, dep, cell in P.cells(GRID, depths=DEPTHS):
        if P.chemistry(sec) != "probe":
            continue
        f = P.scores(cell)
        f["section"], f["depth"] = sec, dep
        rows.append(f)
    S = pd.concat(rows, ignore_index=True)
    S["kind"] = S["kind"].astype(str)
    return S

def cell_values(S, arm, kind, dm, dep):
    if kind == "clean":
        s = S[(S.arm == arm) & (S.condition == "clean") & (S.depth == dep)]
        return s.set_index("section").discard_frac
    s = S[(S.arm == arm) & (S.kind == kind) & (S.diameter_um == dm) & (S.depth == dep)]
    return s.set_index("section").recall

def table(S):
    head = "| | " + " | ".join(c for c, _ in COLS) + " | Sections |"
    out = [head, "|---" * (len(COLS) + 2) + "|"]
    for lab, kind, dm in ROWS:
        cells, ns = [], []
        for dep in DEPTHS:
            run = set(S.section[S.depth == dep])
            v = {c: cell_values(S, f(dm, dep), kind, dm, dep) for c, f in COLS}
            have = [c for c in v if run <= set(v[c].dropna().index)]
            common = sorted(run)
            cells.append({c: (100 * v[c][common].median() if c in have else None) for c in v})
            ns.append(len(common))
        txt = []
        for c, _ in COLS:
            a, b = cells[0][c], cells[1][c]
            txt.append("not yet run" if a is None and b is None else
                       " / ".join("–" if x is None else f"{x:.1f}" for x in (a, b)))
        out.append(f"| {lab} | " + " | ".join(txt) + f" | {ns[0]} / {ns[1]} |")
    fpr, ns = {c: [] for c, _ in COLS}, []
    for dep in DEPTHS:
        common = None
        per = {}
        for c, f in COLS:
            arms = sorted({f(dm, dep) for dm in (32, 64, 256)})
            per[c] = {a: cell_values(S, a, "clean", None, dep).dropna() for a in arms}
            idx = set.intersection(*(set(x.index) for x in per[c].values()))
            common = idx if common is None else common & idx
        common = sorted(common)
        ns.append(len(common))
        for c in per:
            m = [100 * x[common].median() for x in per[c].values()]
            fpr[c].append(f"{min(m):.1f}" if len(m) == 1 else f"{min(m):.1f}–{max(m):.1f}")
    out.append("| Uninjected tissue removed | " + " | ".join(" / ".join(fpr[c]) for c, _ in COLS) +
               f" | {ns[0]} / {ns[1]} |")
    out.append("| Setting depends on the damage diameter | no | yes | yes | no | no | |")
    return "\n".join(out)

def verify():
    sample = [("25_spotsweeper_16um_k48", "dropout01_32um", 100), ("44_spotsweeper_8um_k192", "mito_32um", 12),
              ("19_spatialartifacts_8um", "dropout01_256um", 12), ("20_spatialartifacts_16um", "dropout01edge_64um", 12),
              ("46_spotsweeper_8um_k3072", "mito_256um", 100), ("09_mad_32um", "clean", 12),
              ("19_spatialartifacts_8um", "clean", 100)]
    S = load()
    for arm, cond, dep in sample:
        raw, ref = [], []
        for sec, d, cell in P.cells(GRID, depths=(dep,)):
            if P.chemistry(sec) != "probe":
                continue
            f = os.path.join(cell, "results", arm, f"{cond}.txt")
            if not os.path.exists(f):
                continue
            pred = {ln.strip() for ln in open(f) if ln.strip()}
            d0 = os.path.join(cell, cond)
            gt = pd.read_parquet(os.path.join(d0, "ground_truth_008um.parquet"), columns=["barcode", "is_artifact"])
            uni = set(gt.barcode)
            pos = os.path.join(d0, "square_008um", "spatial", "tissue_positions.parquet")
            if os.path.exists(pos):
                p = pd.read_parquet(pos)
                if "in_tissue" in p.columns:
                    uni = set(p.barcode[p.in_tissue == 1])
            pred &= uni
            if cond == "clean":
                raw.append(len(pred) / len(uni))
            else:
                t = set(gt.barcode[gt.is_artifact]) & uni
                raw.append(len(pred & t) / len(t))
            s = S[(S.section == sec) & (S.depth == dep) & (S.arm == arm) & (S.condition == cond)]
            ref.append(float(s.discard_frac.iloc[0] if cond == "clean" else s.recall.iloc[0]) if len(s) else np.nan)
        diff = np.nanmax(np.abs(np.array(raw) - np.array(ref))) if raw else np.nan
        print(f"{arm:26} {cond:20} d{dep:<3} n={len(raw):2} from barcodes {100*np.median(raw):5.1f}  "
              f"scores.csv {100*np.nanmedian(ref):5.1f}  largest per-section difference {100*diff:.4f} points")

if __name__ == "__main__":
    if "--verify" in sys.argv:
        verify()
    else:
        S = load()
        md = table(S)
        open(os.path.join(GRID, "figures", "guideline_test_table.md"), "w").write(md + "\n")
        print(md)
        for arm in ("19_spatialartifacts_8um", "20_spatialartifacts_16um"):
            for dep in DEPTHS:
                v = cell_values(S, arm, "clean", None, dep).dropna()
                print(f"{arm} d{dep}: removed nothing on the uninjected condition in {int((v == 0).sum())} of {len(v)} sections")

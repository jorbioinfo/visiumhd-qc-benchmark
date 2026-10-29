import argparse, csv, glob, os
import numpy as np, pandas as pd

import panel as P

ARMS = {
    "02_mad_8um": ("section-wide", 8), "08_mad_16um": ("section-wide", 16),
    "09_mad_32um": ("section-wide", 32), "10_mad_64um": ("section-wide", 64),
    **{arm: ("neighbourhood", b) for b, arm in P.SPOTSWEEPER.items()},
    "30_miqc_8um": ("mixture", 8), "31_miqc_64um": ("mixture", 64),
    "32_fixedumi_8um": ("fixed", 8), "33_fixedumi_64um": ("fixed", 64),
    "19_spatialartifacts_8um": ("edge-detector", 8),
    "20_spatialartifacts_16um": ("edge-detector", 16),
}

def split_kind(kind):
    return (kind[:-4], "edge") if kind.endswith("edge") else (kind, "interior")

def cells(grid):
    for sec, depth, path in P.cells(grid):
        yield os.path.basename(path), sec, depth, os.path.join(path, "scored", "scores.csv")

def edge_barcodes(grid, cell, radius_bins):
    from scipy.ndimage import distance_transform_edt
    f = os.path.join(grid, cell, "clean", "square_008um", "spatial",
                     "tissue_positions.parquet")
    if not os.path.exists(f):
        return None
    pos = pd.read_parquet(f)
    tis = pos.in_tissue.to_numpy() == 1
    g = np.zeros((int(pos.array_row.max()) + 1, int(pos.array_col.max()) + 1), bool)
    g[pos.array_row.to_numpy()[tis], pos.array_col.to_numpy()[tis]] = True
    d = distance_transform_edt(g)
    near = d[pos.array_row.to_numpy(), pos.array_col.to_numpy()]
    return set(pos.barcode.to_numpy()[tis & (near > 0) & (near <= radius_bins)]), tis.sum()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--outdir", default="/work")
    ap.add_argument("--edge-radius-bins", type=float, default=16.0,
                    help="half of the largest lesion diameter, in 8 um bins")
    a = ap.parse_args()

    rec, dis, edg = [], [], []
    for cell, ds, depth, p in cells(a.grid):
        df = P.scores(os.path.dirname(os.path.dirname(p)))
        for _, r in df.iterrows():
            if r["arm"] not in ARMS:
                continue
            tool, binum = ARMS[r["arm"]]
            base = dict(dataset=ds, chem=P.chemistry(ds), depth=depth, tool=tool,
                        bin_um=binum)
            if r["kind"] in ("clean", "none") or pd.isna(r.get("diameter_um")):
                if pd.notna(r.get("discard_frac")) and pd.notna(r.get("n_total")):
                    n = int(r["n_total"]); k = int(round(float(r["discard_frac"]) * n))
                    dis.append(dict(base, disc=k, keep=n - k))
                continue
            damage, location = split_kind(str(r["kind"]))
            if pd.isna(r.get("tp")) or (int(r["tp"]) + int(r["fn"])) == 0:
                continue
            rec.append(dict(base, damage=damage, location=location,
                            diam=float(r["diameter_um"]),
                            tp=int(r["tp"]), fn=int(r["fn"])))

        eb = edge_barcodes(a.grid, cell, a.edge_radius_bins)
        if eb is None:
            continue
        edge_set, _ = eb
        for arm, (tool, binum) in ARMS.items():
            f = os.path.join(a.grid, cell, "results", arm, "clean.txt")
            if not os.path.exists(f):
                continue
            got = {ln.strip() for ln in open(f) if ln.strip()}
            k = len(got & edge_set)
            edg.append(dict(dataset=ds, chem=P.chemistry(ds), depth=depth, tool=tool,
                            bin_um=binum,
                            arm=arm, edge_disc=k, edge_keep=len(edge_set) - k,
                            edge_n=len(edge_set)))

    for name, rows in (("recall_table.csv", rec), ("discard_table.csv", dis),
                       ("edge_clean.csv", edg)):
        if not rows:
            print("no rows for", name); continue
        pd.DataFrame(rows).to_csv(os.path.join(a.outdir, name), index=False)
        print(f"wrote {name}: {len(rows)} rows")

if __name__ == "__main__":
    main()

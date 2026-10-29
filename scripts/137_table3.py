import os
import pandas as pd

WK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(WK, "benchmark_grid", "figures")
DEPTHS = (100, 50, 25, 12)
TESTS = (("SpotSweeper", "Neighbourhood (SpotSweeper)"), ("scuttle MAD", "Section-wide (MAD)"))
KINDS = (("dropout01", "signal loss"), ("mito", "mito. elevation"))

d = pd.concat([pd.read_csv(os.path.join(OUT, f"fig5_binsize_d{x}_stats.csv")) for x in DEPTHS])
d = d[d.metric == "recall"]

def table(chem):
    rows = ["| Test | Damage | LQR | 100 | 50 | 25 | 12 UMIs/bin |", "|---|---|---|---:|---:|---:|---:|"]
    gain, ns = 0.0, {}
    for tool, tlab in TESTS:
        for kind, klab in KINDS:
            for i, les in enumerate(sorted(d.lesion_um.unique())):
                cells = []
                for dep in DEPTHS:
                    g = d[(d.chem == chem) & (d.tool == tool) & (d.damage == kind) & (d.lesion_um == les)
                          & (d.depth == dep)].set_index("qc_bin_um").sort_index()
                    best = g.value.max()
                    opt = int(g.index[g.value >= best][0])
                    cells.append(str(opt))
                    ns.setdefault((kind, dep), set()).update(g.n.astype(int))
                    if tool == "scuttle MAD":
                        gain = max(gain, best - g.value.loc[8])
                first = (tlab if kind == KINDS[0][0] and i == 0 else "", klab if i == 0 else "")
                rows.append(f"| {first[0]} | {first[1]} | {int(les)} µm | " + " | ".join(cells) + " |")
    return "\n".join(rows), gain, ns

for chem, name in (("probe", "table3.md"), ("3prime", "tableS5.md")):
    md, gain, ns = table(chem)
    open(os.path.join(OUT, name), "w").write(md + "\n")
    nstr = ", ".join(f"{k} d{dep}: {'/'.join(map(str, sorted(v)))}" for (k, dep), v in sorted(ns.items()))
    print(f"{name}: {chem}; MAD gain at a coarser optimum at most {gain:.1f} points; n {nstr}")

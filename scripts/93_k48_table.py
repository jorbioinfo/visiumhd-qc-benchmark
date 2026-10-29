import csv, collections, os, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import panel as P

DEPTHS = {"d100": "high", "d12": "low"}
KIND = {"dropout01": "signal loss", "mito": "mito. elevation"}
rows = []
for sec, d, path in P.cells("benchmark_grid", (100, 12)):
    for r in P.scores(path).astype(str).to_dict("records"):
        r["cell"], r["dep"] = sec, f"d{d}"
        rows.append(r)

pair = collections.defaultdict(dict)
for r in rows:
    if r["arm"] in (P.SPOTSWEEPER_PACKAGE_DEFAULT_8UM, P.SPOTSWEEPER[8]) and r["kind"] in KIND:
        pair[(r["cell"], r["dep"], r["kind"], float(r["diameter_um"]))][r["arm"]] = r

agg = collections.defaultdict(list)
for key, v in pair.items():
    if len(v) < 2:
        continue
    sec, dep, kind, dia = key
    a = float(v[P.SPOTSWEEPER_PACKAGE_DEFAULT_8UM]["recall"] or 0) * 100
    b = float(v[P.SPOTSWEEPER[8]]["recall"] or 0) * 100
    agg[(P.chemistry(sec), kind, dia, dep)].append((a, b))

print("| Chemistry | Damage | LQR | Depth | *k* = 48, used here | *k* = 36, package default | Difference | Sections |")
print("|---|---|---:|---|---:|---:|---:|---:|")
for (ch, kind, dia, dep) in sorted(agg, key=lambda t: (P.CHEMS.index(t[0]), t[1], t[2], t[3] != "d100")):
    v = agg[(ch, kind, dia, dep)]
    m36 = statistics.median(x[0] for x in v)
    m48 = statistics.median(x[1] for x in v)
    print(f"| {P.CHEM_LABEL[ch]} | {KIND[kind]} | {dia:.0f} µm | {DEPTHS[dep]} | {m48:.1f} | {m36:.1f} | "
          f"{m36 - m48:+.1f} | {len(v)} |")

n = {len(v) for v in agg.values()}
print(f"\nMedian over sections. Sections per row: {sorted(n)}.", file=sys.stderr)

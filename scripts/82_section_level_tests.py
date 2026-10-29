import argparse, csv, math
from collections import defaultdict

import panel as P

BINS = (8, 16, 32, 64)

def sign_p(pos, n):
    if n == 0:
        return float("nan")
    tail = sum(math.comb(n, i) for i in range(pos, n + 1)) / 2 ** n
    return min(1.0, 2 * tail)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--recall", default="/work/recall_table.csv")
    ap.add_argument("--discard", default="/work/discard_table.csv")
    a = ap.parse_args()

    every = [r for r in csv.DictReader(open(a.recall)) if r["location"] == "interior"]
    every_d = list(csv.DictReader(open(a.discard)))
    for ch in P.CHEMS:
        print("#" * 66)
        print(f"# {P.CHEM_LABEL[ch]} sections")
        tests([r for r in every if P.chemistry(r["dataset"]) == ch],
              [r for r in every_d if P.chemistry(r["dataset"]) == ch])

def tests(rows, drows):

    acc = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in rows:
        k = (r["dataset"], r["tool"], float(r["diam"]))
        b = int(r["bin_um"])
        acc[k][b][0] += int(r["tp"]); acc[k][b][1] += int(r["fn"])
    best = {k: max(v, key=lambda b: v[b][0] / max(v[b][0] + v[b][1], 1))
            for k, v in acc.items()}

    for tool in ("neighbourhood", "section-wide"):
        secs = sorted({k[0] for k in best if k[1] == tool})
        print("=" * 66)
        print(f"{tool} — best QC bin size per section (depths and damage types pooled)")
        print(f"{'section':28}" + "".join(f"{str(int(d)) + ' µm':>9}" for d in (32, 64, 256)))
        ge = lt = 0
        for s in secs:
            cells = [best.get((s, tool, float(d))) for d in (32, 64, 256)]
            print(f"{s:28}" + "".join(f"{c if c else '-':>9}" for c in cells))
            if cells[0] and cells[2]:
                ge += cells[2] > cells[0]; lt += cells[2] < cells[0]
        print(f"  coarser optimum for 256 µm than for 32 µm in {ge}/{ge + lt} sections; "
              f"exact two-sided sign test p = {sign_p(ge, ge + lt):.4f}")
        at_ceiling = sum(1 for s in secs if best.get((s, tool, 256.0)) == max(BINS))
        print(f"  at the 64 µm ceiling for 256 µm lesions: {at_ceiling}/{len(secs)} sections")

    d = defaultdict(dict)
    for r in drows:
        d[(r["tool"], int(r["bin_um"]), r["dataset"])][int(r["depth"])] = (
            int(r["disc"]) / (int(r["disc"]) + int(r["keep"])))
    print("=" * 66)
    print("discard rate on undamaged tissue: does it rise as depth falls?")
    for tool, binum in (("section-wide", 8), ("section-wide", 32), ("section-wide", 64),
                        ("neighbourhood", 8), ("neighbourhood", 64)):
        up = dn = 0
        for (t, b, _), v in d.items():
            if (t, b) != (tool, binum):
                continue
            ks = sorted(v)
            if len(ks) < 2:
                continue
            up += v[ks[0]] > v[ks[-1]]; dn += v[ks[0]] < v[ks[-1]]
        print(f"  {tool:14} {binum:2} µm: rises in {up}/{up + dn} sections, "
              f"p = {sign_p(up, up + dn):.4f}")

if __name__ == "__main__":
    main()

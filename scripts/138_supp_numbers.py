import glob, importlib.util, json, math, os
import numpy as np, pandas as pd

WK = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(WK, "benchmark_grid", "figures")
L = []
say = lambda s="": L.append(s)
pct = lambda x: f"{100 * x:.1f}%"
panel = pd.read_csv(os.path.join(OUT, "table1_panel.csv")).set_index("section")
chem = panel.chemistry.to_dict()
_spec = importlib.util.spec_from_file_location("t1", os.path.join(WK, "scripts", "127_table1.py"))
T1 = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(T1)

def sign_p(k, n):
    if n == 0:
        return float("nan")
    k = max(k, n - k)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n)

mk = []
for f in glob.glob(os.path.join(WK, "benchmark_grid", "Visium_HD_*_d*", "scored", "marker.csv")):
    m = pd.read_csv(f)
    run = f.split(os.sep)[-3]
    name, _, dep = run.rpartition("_d")
    m["sec"], m["depth"] = name.replace("Visium_HD_", ""), int(dep)
    mk.append(m)
M = pd.concat(mk)
M = M[M.depth.isin([100, 50, 25, 12])]
M["chem"] = M.sec.map(chem)
M["ratio_naive"] = M.marker_frac_flagged / M.marker_frac_unflagged
say("## S1. Marker check")
ps = M.groupby("sec").panel_share.median().sort_values()
say(f"- panel share of all counts, median over runs: {pct(ps.iloc[0])} ({ps.index[0]}) to {pct(ps.iloc[-1])} "
    f"({ps.index[-1]}), n = {len(ps)} sections")
MAIN = {"02_mad_8um": "MAD 8 um", "10_mad_64um": "MAD 64 um", "24_spotsweeper_8um_k48": "SpotSweeper 8 um",
        "27_spotsweeper_64um_k48": "SpotSweeper 64 um", "19_spatialartifacts_8um": "SpatialArtifacts 8 um",
        "20_spatialartifacts_16um": "SpatialArtifacts 16 um"}
m6 = M[M.arm.isin(MAIN)]
lv = m6.match_scale_um.value_counts(dropna=False)
say(f"- matching level reached, six main arms, all runs and both chemistries: "
    + ", ".join(f"{('not reported' if pd.isna(k) else f'{int(k)} um')}: {v}" for k, v in lv.items())
    + f" of {len(m6)} arm-by-run estimates")
say()
say("## S2. What the removed tissue was (probe-based)")
p = m6[m6.chem == "probe"]
for arm, nm in MAIN.items():
    s = p[p.arm == arm]
    r = s.ratio_matched.dropna()
    d100 = s[s.depth == 100]
    say(f"- {nm}: matched ratio median {r.median():.2f}, below 1 in {int((r < 1).sum())}/{len(r)} runs, range "
        f"{r.min():.2f} to {r.max():.2f}; at 100 UMIs/bin unmatched {d100.ratio_naive.median():.2f} -> matched "
        f"{d100.ratio_matched.median():.2f} (n = {d100.ratio_matched.notna().sum()}); median share removed "
        f"{s.pct_flagged.median():.1f}%; sign p = {sign_p(int((r < 1).sum()), int((r != 1).sum())):.2g}")
t = p[p.arm == "02_mad_8um"].groupby("sec").ratio_matched.median().dropna().sort_values()
say(f"- MAD 8 um matched ratio by section (median over depths): {t.iloc[0]:.2f} ({t.index[0]}) to "
    f"{t.iloc[-1]:.2f} ({t.index[-1]}), n = {len(t)}")
say()

G0 = json.load(open(os.path.join(OUT, "guideline_numbers.json")))

R = pd.read_csv(os.path.join(WK, "recall_table.csv"))
R["recall"] = R.tp / (R.tp + R.fn)
E = pd.read_csv(os.path.join(WK, "edge_clean.csv"))
D = pd.read_csv(os.path.join(WK, "discard_table.csv"))
say("## S3. The tissue boundary (probe-based)")
for dep in (100, 12):
    for dm in (32, 64, 256):
        opts = {}
        for loc in ("interior", "edge"):
            s = R[(R.chem == "probe") & (R.tool == "neighbourhood") & (R.depth == dep) & (R.diam == dm)
                  & (R.location == loc) & R.damage.isin(["dropout01", "mito"])]
            w = s.pivot_table(index=["dataset", "damage"], columns="bin_um", values="recall")
            w = w.groupby(level="dataset").mean().dropna()
            med = w.median()
            opts[loc] = (int(med.idxmax()), len(w))
        say(f"- neighbourhood optimum, {dm} um LQRs d{dep} (mean of the two kinds): interior {opts['interior'][0]} um "
            f"(n = {opts['interior'][1]}), edge {opts['edge'][0]} um (n = {opts['edge'][1]})")
for tool, b in (("section-wide", 8), ("section-wide", 64), ("neighbourhood", 8), ("neighbourhood", 64),
                ("fixed", 8), ("edge-detector", 8), ("edge-detector", 16)):
    for dep in (100, 12):
        e = E[(E.chem == "probe") & (E.tool == tool) & (E.bin_um == b) & (E.depth == dep)].set_index("dataset")
        d = D[(D.chem == "probe") & (D.tool == tool) & (D.bin_um == b) & (D.depth == dep)].set_index("dataset")
        com = e.index.intersection(d.index)
        if not len(com):
            continue
        band = (e.edge_disc / e.edge_n)[com]
        whole = (d.disc / (d.disc + d.keep))[com]
        say(f"- {tool} {b} um d{dep}: undamaged boundary band removed {pct(band.median())}, whole section "
            f"{pct(whole.median())} (n = {len(com)})")
for tool, b in (("section-wide", 8), ("section-wide", 64), ("neighbourhood", 8), ("neighbourhood", 64),
                ("fixed", 8), ("edge-detector", 8), ("edge-detector", 16)):
    for dep in (100, 50, 25, 12):
        out = []
        for dmg in ("dropout01", "mito"):
            for loc in ("interior", "edge"):
                s = R[(R.chem == "probe") & (R.tool == tool) & (R.bin_um == b) & (R.depth == dep) & (R.diam == 256)
                      & (R.damage == dmg) & (R.location == loc)]
                if len(s):
                    out.append(f"{dmg} {loc} {pct(s.recall.median())}")
        if out:
            say(f"  - recall of 256 um LQRs, {tool} {b} um d{dep}: " + ", ".join(out) + f" (n = {len(s)})")
say()

say("## S5. Mitochondrial baseline (per-bin median and 90th percentile, all sections)")
B = pd.read_csv(os.path.join(OUT, "figS4_mito_baseline.csv"))
B["sec"] = B.dataset.str.replace(" ", "_")
B["species"] = B.sec.map(lambda s: panel.species.get(s, "?"))
for ch in ("probe", "3prime"):
    for sp in ("human", "mouse"):
        s = B[(B.chem == ch) & (B.species == sp)].sort_values("median", ascending=False)
        if len(s):
            say(f"- {ch} {sp} (n = {len(s)}): highest " + ", ".join(f"{r.dataset} {r['median']:.1f}% (p90 {r.p90:.1f}%)"
                for _, r in s.head(3).iterrows()) + f"; lowest {s.iloc[-1].dataset} {s.iloc[-1]['median']:.1f}%")
pr = B[B.sec == "Human_Prostate_Cancer_FFPE"]
if len(pr):
    say(f"- human prostate cancer FFPE: per-bin median {pr['median'].iloc[0]:.1f}%, native mean UMIs/bin "
        f"{panel.native_mean_umi.get('Human_Prostate_Cancer_FFPE', float('nan')):.0f}, rank by native depth "
        f"{int(panel.native_mean_umi.rank().get('Human_Prostate_Cancer_FFPE'))} of {len(panel)}")
say(f"- section mean mitochondrial baseline (Table 1): human kidney FFPE "
    f"{panel.mito_baseline_pct.get('Human_Kidney_FFPE'):.2f}%, panel median {panel.mito_baseline_pct.median():.2f}%")
for tool, b in (("section-wide", 8), ("neighbourhood", 8), ("neighbourhood", 64)):
    for dm in (32, 64, 256):
        s = R[(R.tool == tool) & (R.bin_um == b) & (R.depth == 100) & (R.diam == dm) & (R.damage == "mito")
              & (R.location == "interior") & (R.chem == "probe")].set_index("dataset").recall
        if "Human_Kidney_FFPE" in s.index:
            say(f"  - mito {dm} um d100, {tool} {b} um: human kidney FFPE {pct(s['Human_Kidney_FFPE'])}, "
                f"other probe sections median {pct(s.drop('Human_Kidney_FFPE').median())} (n = {len(s) - 1})")
for dm in (32, 64, 256):
    s = R[(R.tool == "section-wide") & (R.bin_um == 8) & (R.depth == 100) & (R.diam == dm) & (R.damage == "dropout01")
          & (R.location == "interior")].set_index("dataset").recall
    say(f"  - signal loss {dm} um d100, section-wide 8 um: human kidney FFPE {pct(s['Human_Kidney_FFPE'])}")
say()

St = json.load(open(os.path.join(OUT, "guideline_settings.json")))
NB, VL, WC = G0["neighbourhood"], G0["very_large"], G0["widened_check"]
LQR = {0: 32, 1: 64, 2: 256, 3: 512}
rows = ["| Damage size | UMIs/bin | Setting | Signal loss | Mito. elevation | Sections |",
        "|---|---|---|---:|---:|---:|"]
for i, c in enumerate(St["classes"]):
    dm = LQR[i]
    for tier, dep in (("hi", 100), ("lo", 12)):
        st = c[tier]
        if st["bin"] == 8 and st["k"] > 48:
            e = WC[f"{dm}|{dep}"]; key = f"k{st['k']}"
            vals = [e[k]["median"][key] for k in ("dropout01", "mito")]; n = e["dropout01"]["n"]
        elif dm == 512:
            vals = [VL[f"{k}|{dep}"]["median"][str(st["bin"])] for k in ("dropout01", "mito")]
            n = VL[f"dropout01|{dep}"]["n"]
        else:
            vals = [NB[f"probe|{k}|{dm}|{dep}"]["median"][str(st["bin"])] for k in ("dropout01", "mito")]
            n = NB[f"probe|dropout01|{dm}|{dep}"]["n"]
        setting = f"{st['bin']} µm bins" + ("" if st["bin"] == 8 else " + A&D") + f", {st['k']} neighbours"
        size = f"{c['label']} ({dm} µm LQRs)" if tier == "hi" else ""
        rows.append(f"| {size} | {'≥ 100' if tier == 'hi' else '12'} | "
                    f"{setting} | {100 * vals[0]:.1f}% | {100 * vals[1]:.1f}% | {n} |")
open(os.path.join(OUT, "guideline_table.md"), "w").write("\n".join(rows) + "\n")
say("## S6: guideline_table.md written")
say()

say("## S9. 36 against 48 neighbours, cost on undamaged tissue (probe)")
for dep in (100, 12):
    a = D[(D.chem == "probe") & (D.tool == "neighbourhood") & (D.bin_um == 8) & (D.depth == dep)].set_index("dataset")
    a = a.disc / (a.disc + a.keep)
    b = []
    for f in glob.glob(os.path.join(WK, "benchmark_grid", f"Visium_HD_*_d{dep}", "scored", "scores.csv")):
        s = pd.read_csv(f)
        s = s[(s.arm == "03_spotsweeper_8um") & (s.condition == "clean")]
        if len(s):
            b.append((f.split(os.sep)[-3].replace("Visium_HD_", "").rpartition("_d")[0], s.discard_frac.iloc[0]))
    b = pd.Series(dict(b))
    com = a.index.intersection(b.index).intersection([s for s, c in chem.items() if c == "probe"])
    say(f"- d{dep}: 48 neighbours {pct(a[com].median())}, 36 neighbours {pct(b[com].median())} (n = {len(com)})")
say()

PAIRS = [
    ("Mouse brain", ["Mouse_Brain", "Mouse_Brain_Fresh_Frozen", "Mouse_Brain_Fixed_Frozen", "3prime_Mouse_Brain"]),
    ("Human tonsil", ["Human_Tonsil_Fresh_Frozen", "Human_Tonsil_Fresh_Frozen_IF", "3prime_Human_Tonsil"]),
    ("Human lung cancer", ["Human_Lung_Cancer", "Human_Lung_Cancer_Fixed_Frozen", "SixCancer_Lung_Cancer"]),
    ("Human lymph node", ["Human_Lymph_Node_FFPE", "3prime_Human_Lymph_Node"]),
    ("Human ovarian cancer", ["Human_Ovarian_Cancer_FF", "3prime_Human_Ovarian_Cancer", "SixCancer_Ovarian_Cancer"]),
    ("Mouse kidney", ["Mouse_Kidney", "3prime_Mouse_Kidney"]),
    ("Mouse embryo", ["Mouse_Embryo", "11mm_Mouse_Embryo", "3prime_Mouse_Embryo"]),
    ("Human pancreas", ["Human_Pancreas", "11mm_Human_Pancreas"]),
    ("Human colon cancer", ["Human_Colon_Cancer", "11mm_Human_Colon_Cancer_Fresh_Frozen", "Human_Colon_Cancer_P1",
                            "Human_Colon_Cancer_P5", "SixCancer_Colorectal_Cancer"]),
    ("Human breast cancer", ["Human_Breast_Cancer_Fresh_Frozen", "11mm_Human_Breast_Cancer", "Human_Breast_Cancer_FFPE"]),
]
nb = R[(R.tool == "neighbourhood") & (R.location == "interior") & R.damage.isin(["dropout01", "mito"])]

def optimum(sec, dm):
    s = nb[(nb.dataset == sec) & (nb.diam == dm)]
    if s.empty:
        return "–"
    return str(int(s.groupby("bin_um").recall.mean().idxmax()))

mad8 = D[(D.tool == "section-wide") & (D.bin_um == 8)].copy()
mad8["cost"] = mad8.disc / (mad8.disc + mad8.keep)
rows = ["| Tissue | Section | Chemistry, preparation | Space Ranger | Optimum bin for 32 / 64 / 256 µm LQRs (µm) "
        "| MAD cost on 8 µm bins, 100 → 12 UMIs/bin |", "|---|---|---|---|---|---|"]
for tissue, secs in PAIRS:
    for i, sec in enumerate(secs):
        c = mad8[mad8.dataset == sec].set_index("depth").cost
        cost = (f"{pct(c[100])} → {pct(c[12])}" if 100 in c.index and 12 in c.index
                else f"no 100 level; {pct(c[50])} at 50 → {pct(c[12])}" if 50 in c.index else "–")
        ch = "3′ polyA" if chem.get(sec) == "3prime" else "probe"
        r1 = {"section": sec, "chemistry": chem.get(sec)}
        sr = panel.spaceranger.get(sec)
        sr = sr if isinstance(sr, str) and sr else T1.SPACERANGER.get(sec, "-")
        slide = "11 mm" if panel.slide_mm.get(sec) == 11.0 else "6.5 mm"
        rows.append(f"| {tissue if i == 0 else ''} | {T1.label(sec)}{', ' + slide if slide == '11 mm' else ''} | "
                    f"{ch}, {T1.prep(r1)} | {sr} | "
                    f"{' / '.join(optimum(sec, dm) for dm in (32, 64, 256))} | {cost} |")
open(os.path.join(OUT, "tableS4.md"), "w").write("\n".join(rows) + "\n")
say(f"## S10: Table S4 written, {sum(len(s) for _, s in PAIRS)} sections in {len(PAIRS)} tissues")

G = json.load(open(os.path.join(OUT, "guideline_numbers.json")))
cols = ["b8", "b16", "b32", "b64", "k192", "k768", "k3072"]
rows = ["| LQR | UMIs/bin | Damage | n | 8 µm | 16 µm | 32 µm | 64 µm | 192 nb. | 768 nb. | 3072 nb. |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
for key, e in G["widened_check"].items():
    dm, dep = key.split("|")
    for k, lab in (("dropout01", "signal loss"), ("mito", "mito. elevation")):
        m = e[k]["median"]
        rows.append(f"| {dm} µm | {dep} | {lab} | {e[k]['n']} | " + " | ".join(f"{100 * m[c]:.1f}" for c in cols) + " |")
open(os.path.join(OUT, "tableS1.md"), "w").write("\n".join(rows) + "\n")
say(f"## Table S1 written, {len(rows) - 2} rows")

src = os.path.join(WK, "data", "sources.tsv")
if os.path.exists(src):
    S = pd.read_csv(src, sep="\t")
    rows = ["| Section | Provider | Source | Licence |", "|---|---|---|---|"]
    for _, r in S.iterrows():
        lic = str(r.licence).split(" (")[0]
        rows.append(f"| {r.section.replace('Visium_HD_', '').replace('_', ' ')} | {r.provider} | {r.source} | {lic} |")
    open(os.path.join(OUT, "tableS6.md"), "w").write("\n".join(rows) + "\n")
    say(f"## Table S6 written, {len(S)} sections; providers {S.provider.value_counts().to_dict()}")

open(os.path.join(OUT, "supp_numbers.md"), "w").write("# Supplementary numbers\n\n" + "\n".join(L) + "\n")
print("wrote supp_numbers.md, tableS4.md, tableS6.md;", len(L), "lines")

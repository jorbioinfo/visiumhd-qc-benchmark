import glob, json, math, os
import numpy as np, pandas as pd
from scipy import stats

import panel as P

GRID = "benchmark_grid"
OUT = os.path.join(GRID, "figures")
SS, MAD = P.SPOTSWEEPER, P.MAD
lines = []

def say(s=""):
    lines.append(s)

def sign_p(k, n):
    if n == 0:
        return float("nan")
    k = max(k, n - k)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k, n + 1)) / 2 ** n)

def pct(x):
    return f"{100 * x:.1f}%"

def med(s):
    return float(np.median(s)) if len(s) else float("nan")

def md_table(df, index=True):
    cols = ([df.index.name or ""] if index else []) + [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for i, r in df.iterrows():
        out.append("| " + " | ".join(([str(i)] if index else []) + [str(v) for v in r]) + " |")
    return "\n".join(out)

def panel_table():
    rows = []
    for sec, (sp, chem, mm) in sorted(P.SECTIONS.items()):
        runs = sorted(glob.glob(f"{GRID}/{P.PREFIX}{sec}_d*"))
        depths = sorted(int(r.rsplit("_d", 1)[1]) for r in runs)
        native = base = None
        for r in runs:
            try:
                ip = json.load(open(f"{r}/injection_params.json"))
                native = ip["target_umi_per_bin"] / ip["subsample"]
                core = f"{r}/injection_params_mito_relative_core.json"
                base = (json.load(open(core)) if os.path.exists(core) else ip).get("mito_baseline_value")
                break
            except (OSError, KeyError, ValueError):
                continue
        n_tissue = None
        for r in runs:
            f = f"{r}/clean/square_008um/spatial/tissue_positions.parquet"
            if os.path.exists(f):
                n_tissue = int((pd.read_parquet(f, columns=["in_tissue"]).in_tissue == 1).sum())
                break
        sr = None
        try:
            import h5py
            h = f"public/{P.PREFIX}{sec}/binned_outputs/square_008um/filtered_feature_bc_matrix.h5"
            with h5py.File(h, "r") as fh:
                v = fh.attrs.get("software_version")
                sr = (v.decode() if isinstance(v, bytes) else str(v)).replace("spaceranger-", "") if v is not None else None
        except Exception:
            pass
        rows.append(dict(section=sec, species=sp, chemistry=chem, slide_mm=mm,
                         depths=",".join(map(str, [d for d in depths if d in P.DEPTHS])),
                         native_mean_umi=round(native, 1) if native else None,
                         mito_baseline_pct=round(100 * base, 2) if base is not None else None,
                         in_tissue_bins=n_tissue, spaceranger=sr))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(OUT, "table1_panel.csv"), index=False)
    return t

t1 = panel_table()
say("# Manuscript numbers")
say()
say(f"Generated {pd.Timestamp.now(tz='Asia/Taipei'):%Y-%m-%d %H:%M} by scripts/126_manuscript_numbers.py.")
say()
say("## Panel")
say(f"- sections: {len(t1)}; probe {sum(t1.chemistry == P.PROBE)}, 3' polyA {sum(t1.chemistry == P.POLYA)}; "
    f"human {sum(t1.species == 'human')}, mouse {sum(t1.species == 'mouse')}; 11 mm slide {sum(t1.slide_mm == 11.0)}")
for ch in P.CHEMS:
    s = t1[t1.chemistry == ch]
    say(f"- native mean UMIs per 8 um bin, {ch}: {s.native_mean_umi.min():.0f} to {s.native_mean_umi.max():.0f} "
        f"(median {s.native_mean_umi.median():.0f}); mito baseline mean {s.mito_baseline_pct.min():.2f}% to "
        f"{s.mito_baseline_pct.max():.2f}% (median {s.mito_baseline_pct.median():.2f}%)")
say(f"- native below 100 UMIs/bin (no 100 level): {', '.join(t1[t1.native_mean_umi < 100].section)}")
cells = list(P.cells(GRID))
by_d = pd.Series([d for _, d, _ in cells]).value_counts().sort_index(ascending=False)
say(f"- scored runs at the four levels: {len(cells)} ({', '.join(f'{d}: {n}' for d, n in by_d.items())})")

rows = []
for sec, dep, cell in cells:
    f = P.scores(cell)
    f["section"], f["depth"], f["chem"] = sec, dep, P.chemistry(sec)
    rows.append(f)
S = pd.concat(rows, ignore_index=True)
S["kind"] = S["kind"].astype(str)
mito_runs = S[S.kind == "mito"].groupby(["section", "depth"]).size()
say(f"- runs with mitochondrial conditions on the relative rule so far: {len(mito_runs)} of {len(cells)}; "
    f"by depth {dict(mito_runs.reset_index().depth.value_counts().sort_index(ascending=False))}")
edge_runs = S[S.kind.str.endswith("edge")].groupby(["section", "depth"]).size()
say(f"- runs with edge conditions so far: {len(edge_runs)} of {len(cells)}")
say()

def recall(arm, kind, diam, depth, chem):
    s = S[(S.arm == arm) & (S.kind == kind) & (S.diameter_um == diam) & (S.depth == depth) & (S.chem == chem)]
    return s.set_index("section").recall.dropna()

def disc(arm, depth, chem):
    s = S[(S.arm == arm) & (S.condition == "clean") & (S.depth == depth) & (S.chem == chem)]
    return s.set_index("section").discard_frac.dropna()

KIND = {"dropout01": "signal loss", "mito": "mito elevation"}

say("## Failure 1: neighbourhood test (SpotSweeper, 48 neighbours)")
for ch in P.CHEMS:
    for dep in (100, 12):
        for k, kl in KIND.items():
            parts = []
            for dm in (32, 64, 256):
                r8 = recall(SS[8], k, dm, dep, ch)
                parts.append(f"{dm} um {pct(med(r8))}" if len(r8) else f"{dm} um -")
            n = len(recall(SS[8], k, 256, dep, ch))
            say(f"- native 8 um recall, {ch} d{dep} {kl} (n = {n}): " + ", ".join(parts))
say()
for ch in P.CHEMS:
    for dep in (100, 12):
        for k, kl in KIND.items():
            for dm in (64, 256):
                a, b = recall(SS[8], k, dm, dep, ch), recall(SS[64], k, dm, dep, ch)
                j = a.index.intersection(b.index)
                if not len(j):
                    continue
                up = int((b[j] > a[j]).sum()); dn = int((b[j] < a[j]).sum())
                say(f"- A&D 64 um vs native, {ch} d{dep} {kl} {dm} um: median {pct(med(a[j]))} -> {pct(med(b[j]))}; "
                    f"higher in {up}/{up + dn} (sign p = {sign_p(up, up + dn):.3g}); n = {len(j)}")
say()
BINS = (8, 16, 32, 64)
say("### Bin size at which recall peaks (median over sections), SpotSweeper and MAD")
opt_rows = []
for tool, arms in (("neighbourhood", SS), ("section-wide", MAD)):
    for ch in P.CHEMS:
        for k, kl in KIND.items():
            for dm in (32, 64, 256):
                cells_ = []
                for dep in P.DEPTHS:
                    curves = {b: recall(arms[b], k, dm, dep, ch) for b in BINS}
                    common = set.intersection(*(set(c.index) for c in curves.values()))
                    if not common:
                        cells_.append("-"); continue
                    m = [med(curves[b][sorted(common)]) for b in BINS]
                    cells_.append(f"{BINS[int(np.nanargmax(m))]} (n={len(common)})")
                opt_rows.append(dict(tool=tool, chem=ch, damage=kl, lqr=dm, **dict(zip(map(str, P.DEPTHS), cells_))))
say(md_table(pd.DataFrame(opt_rows), index=False))
say()
say("### Per section: optimum for 256 um coarser than for 32 um (depths pooled)")
for tool, arms in (("neighbourhood", SS), ("section-wide", MAD)):
    for ch in P.CHEMS:
        for k, kl in KIND.items():
            gt = lt = eq = 0; secs = 0
            vals = {}
            for sec in sorted(S[S.chem == ch].section.unique()):
                o = {}
                for dm in (32, 256):
                    m = []
                    for b in BINS:
                        v = [recall(arms[b], k, dm, dep, ch).get(sec, np.nan) for dep in P.DEPTHS]
                        m.append(np.nanmean(v) if not all(np.isnan(v)) else np.nan)
                    if all(np.isnan(m)):
                        break
                    o[dm] = BINS[int(np.nanargmax(m))]
                if len(o) < 2:
                    continue
                secs += 1; vals[sec] = (o[32], o[256])
                gt += o[256] > o[32]; lt += o[256] < o[32]; eq += o[256] == o[32]
            pairs = pd.Series(vals).value_counts()
            say(f"- {tool}, {ch}, {kl}: coarser in {gt}, finer in {lt}, tied in {eq} of {secs} sections "
                f"(sign p = {sign_p(gt, gt + lt):.3g}); (32 um, 256 um) optimum pairs: "
                + ", ".join(f"{a}/{b}: {n}" for (a, b), n in pairs.items()))
say()

say("## Widened neighbourhood against A&D (probe, sections with every setting)")
cand = {**{f"b{b}": SS[b] for b in (8, 16, 32, 64)}, **{f"k{48 * (b // 8) ** 2}": a for b, a in P.WIDENED.items()}}
for dm in (32, 64, 256):
    for dep in (100, 12):
        per = {}
        for k in KIND:
            vals = {nm: recall(a, k, dm, dep, "probe") for nm, a in cand.items()}
            com = sorted(set.intersection(*(set(v.index) for v in vals.values())))
            per[k] = {nm: med(v[com]) for nm, v in vals.items()}
            say(f"- {dm} um {KIND[k]} d{dep} (n = {len(com)}): " + ", ".join(f"{nm} {pct(x)}" for nm, x in per[k].items()))
        mean = {nm: np.mean([per[k][nm] for k in KIND]) for nm in cand}
        bb = max((nm for nm in mean if nm.startswith("b")), key=mean.get)
        bk = max((nm for nm in mean if nm.startswith("k")), key=mean.get)
        say(f"  - mean of the two kinds, {dm} um d{dep}: " + ", ".join(f"{nm} {pct(x)}" for nm, x in mean.items())
            + f"; best bin {bb} {pct(mean[bb])}, best widened {bk} {pct(mean[bk])}")
say()
say("## Edge LQRs as the text quotes them (probe)")
for dep in (100, 12):
    e = [med(recall(SS[64], k + "edge", 256, dep, "probe")) for k in KIND]
    i = [med(recall(SS[64], k, 256, dep, "probe")) for k in KIND]
    say(f"- SpotSweeper 64 um, 256 um LQRs d{dep}, mean of the two kinds' medians: edge {pct(np.mean(e))} "
        f"(n = {len(recall(SS[64], 'dropout01edge', 256, dep, 'probe'))}), interior {pct(np.mean(i))}")
SA = {8: "19_spatialartifacts_8um", 16: "20_spatialartifacts_16um"}
for b, a in SA.items():
    for dep in (100, 12):
        r = {dm: recall(a, "dropout01edge", dm, dep, "probe") for dm in (64, 256)}
        say(f"- SpatialArtifacts {b} um, edge signal loss d{dep}: 64 um {pct(med(r[64]))}, 256 um {pct(med(r[256]))} "
            f"(n = {len(r[256])})")
mx = max(((med(recall(a, "mitoedge", dm, dep, "probe")), b, dm, dep) for b, a in SA.items()
          for dm in (32, 64, 256) for dep in (100, 12)))
say(f"- SpatialArtifacts, highest median recall of any edge mitochondrial LQR: {pct(mx[0])} "
    f"({mx[1]} um bins, {mx[2]} um LQRs, d{mx[3]})")
say()

say("## Very large LQRs (512 um), SpotSweeper on 64 and 128 um bins with A&D")
SS128 = "28_spotsweeper_128um_k48"
vl_secs = S[S.diameter_um == 512].groupby("chem").section.nunique()
say(f"- sections holding 512 um LQRs: {int(vl_secs.sum())} ({', '.join(f'{c}: {n}' for c, n in vl_secs.items())}); "
    f"runs: {S[S.diameter_um == 512].groupby(['section', 'depth']).ngroups}")
for ch in P.CHEMS:
    for k in KIND:
        for dep in (100, 12):
            a, b = recall(SS[64], k, 512, dep, ch), recall(SS128, k, 512, dep, ch)
            com = a.index.intersection(b.index)
            if not len(com):
                continue
            dlt = b[com] - a[com]
            up, nz = int((dlt > 0).sum()), int((dlt != 0).sum())
            say(f"- 512 um {KIND[k]}, {ch} d{dep}: 64 um {pct(med(a[com]))} -> 128 um {pct(med(b[com]))}; "
                f"higher in {up}/{nz} (ties {len(com) - nz}; sign p = {sign_p(up, nz):.2g}); n = {len(com)}")
for ch in P.CHEMS:
    for dep in (100, 12):
        a, b = disc(SS[64], dep, ch), disc(SS128, dep, ch)
        com = a.index.intersection(b.index)
        if len(com):
            say(f"- cost on the uninjected condition, {ch} d{dep}: SpotSweeper 64 um {pct(med(a[com]))}, "
                f"128 um {pct(med(b[com]))} (n = {len(com)})")
say()

say("## Failure 2: section-wide MAD on undamaged tissue (clean condition)")
for ch in P.CHEMS:
    for b in BINS:
        meds = []
        for dep in P.DEPTHS:
            d = disc(MAD[b], dep, ch)
            meds.append(f"d{dep} {pct(med(d))} (n={len(d)})")
        say(f"- MAD {b} um, {ch}: " + ", ".join(meds))
    for b in BINS:
        a, c = disc(MAD[b], 100, ch), disc(MAD[b], 12, ch)
        j = a.index.intersection(c.index)
        if not len(j):
            continue
        up = int((c[j] > a[j]).sum()); dn = int((c[j] < a[j]).sum())
        ratio = (c[j] / a[j].replace(0, np.nan)).dropna()
        dpp = (c[j] - a[j]) * 100
        say(f"  - MAD {b} um, {ch}, d12 vs d100 within section (n = {len(j)}): rises in {up}/{up + dn} "
            f"(sign p = {sign_p(up, up + dn):.3g}); ratio {ratio.min():.2f} to {ratio.max():.2f} "
            f"(median {ratio.median():.2f}); change {dpp.min():+.1f} to {dpp.max():+.1f} points "
            f"(median {dpp.median():+.1f})")
say()
parts = []
for b in BINS:
    k = n = 0
    for ch in P.CHEMS:
        a, c = disc(MAD[b], 100, ch), disc(MAD[b], 12, ch)
        j = a.index.intersection(c.index)
        n += len(j); k += int((c[j] > 1.25 * a[j]).sum())
    parts.append(f"{b} um {k}/{n}")
say("- MAD, sections whose cost at d12 exceeds 1.25 x cost at d100: " + ", ".join(parts))
say()
for name, arm in (("SpotSweeper 8 um", SS[8]), ("SpotSweeper 64 um", SS[64]),
                  ("SpatialArtifacts 8 um", "19_spatialartifacts_8um"), ("fixed 100 UMIs 8 um", "32_fixedumi_8um")):
    for ch in P.CHEMS:
        meds = [f"d{dep} {pct(med(disc(arm, dep, ch)))}" for dep in P.DEPTHS]
        a, c = disc(arm, 100, ch), disc(arm, 12, ch)
        j = a.index.intersection(c.index)
        up = int((c[j] > a[j]).sum()); dn = int((c[j] < a[j]).sum())
        say(f"- {name}, {ch}: " + ", ".join(meds) + f"; rises d100 to d12 in {up}/{up + dn} (n = {len(j)})")
say()

say("### F1 for 256 um LQRs, MAD, median over sections")
for ch in P.CHEMS:
    for k, kl in KIND.items():
        tab = {}
        for b in BINS:
            row = {}
            for dep in P.DEPTHS:
                s = S[(S.arm == MAD[b]) & (S.kind == k) & (S.diameter_um == 256) & (S.depth == dep) & (S.chem == ch)]
                row[f"d{dep}"] = f"{med(s.f1.dropna()):.3f} (n={s.f1.notna().sum()})" if len(s) else "-"
            tab[f"{b} um"] = row
        say(f"{ch}, {kl}:")
        say(md_table(pd.DataFrame(tab).T))
        say()

say("### miQC")
for arm in ("30_miqc_8um", "31_miqc_64um"):
    fails = P.method_failures(GRID, arm, P.big512_runs("."))
    if any(f for f, _ in fails.values()):
        say(f"- {arm} failed verdicts by depth (failed/owed): " +
            ", ".join(f"d{d} {f}/{n} ({pct(f / n) if n else '-'})" for d, (f, n) in fails.items()))
for ch in P.CHEMS:
    for k, kl in KIND.items():
        for arm in ("31_miqc_64um", MAD[64]):
            s = S[(S.arm == arm) & (S.kind == k) & (S.diameter_um == 256) & (S.chem == ch) & (S.depth.isin([100, 12]))]
            if len(s):
                say(f"- {arm}, {ch}, {kl} 256 um, F1 at d100 / d12: "
                    f"{med(s[s.depth == 100].f1.dropna()):.3f} / {med(s[s.depth == 12].f1.dropna()):.3f} "
                    f"(n = {s[s.depth == 100].f1.notna().sum()} / {s[s.depth == 12].f1.notna().sum()})")
for ch in P.CHEMS:
    for dep in (100, 12):
        say(f"- miQC 64 um discard on clean, {ch} d{dep}: {pct(med(disc('31_miqc_64um', dep, ch)))} "
            f"(n = {len(disc('31_miqc_64um', dep, ch))}); MAD 64 um {pct(med(disc(MAD[64], dep, ch)))}")
say()

say("### Tissue boundary (uninjected boundary band against the whole section)")
try:
    E = pd.read_csv("edge_clean.csv"); D = pd.read_csv("discard_table.csv")
    E["edge_rate"] = E.edge_disc / E.edge_n
    D["rate"] = D.disc / (D.disc + D.keep)
    M = E.merge(D, on=["dataset", "chem", "depth", "tool", "bin_um"])
    for ch in P.CHEMS:
        for b in BINS:
            m = M[(M.tool == "section-wide") & (M.bin_um == b) & (M.chem == ch) & (M.rate > 0)]
            r = m.edge_rate / m.rate
            say(f"- MAD {b} um, {ch}: boundary / whole-section discard ratio median {r.median():.2f} "
                f"(n = {len(m)} runs)")
except OSError:
    say("- edge_clean.csv not found")
say()

say("### Edge LQRs (subset so far)")
for ch in P.CHEMS:
    for k, kl in (("dropout01edge", "signal loss"), ("mitoedge", "mito elevation")):
        for arm, nm in ((SS[8], "SpotSweeper 8"), (SS[64], "SpotSweeper 64"), (MAD[8], "MAD 8"), (MAD[64], "MAD 64")):
            s = S[(S.arm == arm) & (S.kind == k) & (S.diameter_um == 256) & (S.depth == 100) & (S.chem == ch)]
            if len(s):
                say(f"- {nm}, {ch}, edge {kl} 256 um d100: recall median {pct(med(s.recall.dropna()))} (n = {len(s)})")
say()

say("## Downstream clustering (32 um, Leiden), clean condition")
dsb = pd.read_csv(os.path.join(OUT, "downstream_binned.csv")).drop_duplicates(["section", "arm", "comparison"], keep="last")
dsb["sec"] = dsb.section.map(P.short); dsb["chem"] = dsb.sec.map(P.chemistry)
rs = dsb[dsb.comparison == "reseed"]
say(f"- sections: {dsb.sec.nunique()} ({', '.join(f'{c}: {dsb[dsb.chem == c].sec.nunique()}' for c in P.CHEMS)})")
say(f"- reseed ARI (no QC, same depth): {rs.ari.min():.3f} to {rs.ari.max():.3f}, median {rs.ari.median():.3f}")
ad = dsb[dsb.comparison == "across_depth"].pivot(index="sec", columns="arm", values="ari")
say(f"- no-QC ARI across depth: median {ad['01_no_qc'].median():.3f} ({ad['01_no_qc'].min():.3f} to {ad['01_no_qc'].max():.3f})")
for ch in P.CHEMS:
    sub = ad[ad.index.map(P.chemistry) == ch]
    for arm in [c for c in ad.columns if c != "01_no_qc"]:
        dlt = (sub[arm] - sub["01_no_qc"]).dropna()
        lo = int((dlt < 0).sum()); hi = int((dlt > 0).sum())
        w = stats.wilcoxon(dlt).pvalue if len(dlt) > 5 else float("nan")
        say(f"- {arm}, {ch} (n = {len(dlt)}): lower than no QC in {lo}/{lo + hi}, median {dlt.median():+.3f} "
            f"({dlt.min():+.3f} to {dlt.max():+.3f}); sign p = {sign_p(lo, lo + hi):.2g}, Wilcoxon p = {w:.2g}")
swing = []
for sec in ad.index:
    a, c = disc(MAD[8], 100, P.chemistry(sec)).get(sec), disc(MAD[8], 12, P.chemistry(sec)).get(sec)
    if a is not None and c is not None and "02_mad_8um" in ad.columns:
        swing.append(dict(sec=sec, chem=P.chemistry(sec), swing=100 * (c - a),
                          d=ad.loc[sec, "02_mad_8um"] - ad.loc[sec, "01_no_qc"]))
sw = pd.DataFrame(swing).dropna()
for ch in P.CHEMS:
    s = sw[sw.chem == ch]
    if len(s) > 3:
        rho, p = stats.spearmanr(s.swing, s.d)
        say(f"- native MAD discard swing (d12 - d100, points) against ARI change, {ch} (n = {len(s)}): "
            f"Spearman rho {rho:.2f}, p = {p:.2g}; swing {s.swing.min():.1f} to {s.swing.max():.1f}")
rsd = rs.set_index(["sec", "arm"]).ari
for sec in ad.index[(ad.index.map(P.chemistry) == "probe")]:
    if ad.loc[sec, "10_mad_64um"] < ad.loc[sec, "02_mad_8um"] - 0.1:
        c = {d: disc(MAD[64], d, "probe").get(sec, float("nan")) for d in (100, 12)}
        say(f"- {sec}: ARI no QC {ad.loc[sec, '01_no_qc']:.3f}, native MAD {ad.loc[sec, '02_mad_8um']:.3f}, "
            f"MAD 64 um {ad.loc[sec, '10_mad_64um']:.3f}; MAD 64 um cost d100 {pct(c[100])}, d12 {pct(c[12])}; "
            f"reseed ARI of MAD 64 um survivors {rsd.get((sec, '10_mad_64um'), float('nan')):.3f}")
sw.to_csv(os.path.join(OUT, "downstream_swing.csv"), index=False)
say()

say("## Marker check (uninjected, matched on location and library size)")
mk = []
for sec, dep, cell in P.cells(GRID, need=("scored", "marker.csv")):
    m = pd.read_csv(f"{cell}/scored/marker.csv"); m["sec"], m["depth"], m["chem"] = sec, dep, P.chemistry(sec)
    mk.append(m)
if mk:
    MK = pd.concat(mk)
    for arm in (MAD[8], MAD[64], SS[8]):
        for ch in P.CHEMS:
            s = MK[(MK.arm == arm) & (MK.chem == ch)].ratio_matched.dropna()
            if len(s):
                say(f"- {arm}, {ch}: matched marker ratio, discarded / kept, median {s.median():.2f} "
                    f"(below 1 in {int((s < 1).sum())}/{len(s)} runs)")
say()

open(os.path.join(OUT, "manuscript_numbers.md"), "w").write("\n".join(lines) + "\n")
print("wrote", os.path.join(OUT, "manuscript_numbers.md"), len(lines), "lines")

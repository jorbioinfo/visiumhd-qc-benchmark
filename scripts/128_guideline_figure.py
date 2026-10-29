import importlib.util, json, math, os, sys
import numpy as np, pandas as pd

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WK, "scripts"))
import panel as P

GRID = os.path.join(WK, "benchmark_grid")
OUT = os.path.join(GRID, "figures")
SS, MAD, WIDE = P.SPOTSWEEPER, P.MAD, P.WIDENED
BINS = (8, 16, 32, 64)
KINDS = ("dropout01", "mito")
LQRS = (32, 64, 256)
ENDS = (100, 12)
TOL = 1.25
WIDE_MIN = 4
VERY_LARGE = 512
BINS_VL = BINS + (128,)

COVERAGE_ARMS = {
    "01_no_qc": "left out: reference, removes nothing",
    "02_mad_8um": "MAD outcome, at or above the depth cut",
    "08_mad_16um": "MAD bin rule: cost grows as depth falls",
    "09_mad_32um": "MAD outcome, below the depth cut",
    "10_mad_64um": "MAD bin rule: coarser than needed; boundary corrections compared in 129",
    "03_spotsweeper_8um": "Discussion, avoid: the package default of 36 neighbours",
    "04_spotsweeper_16um": "left out: 36-neighbour run of 25, superseded",
    "05_spotsweeper_32um": "left out: 36-neighbour run of 26, superseded",
    "06_spotsweeper_64um": "left out: 36-neighbour run of 27, superseded",
    "19_spatialartifacts_8um": "edge outcome, at or above the depth cut",
    "20_spatialartifacts_16um": "edge outcome, below the depth cut",
    "21_spotsweeper_8um_k144": "left out: widened from 36 neighbours, superseded by 44",
    "22_spotsweeper_8um_k576": "left out: widened from 36 neighbours, superseded by 45",
    "23_spotsweeper_8um_k2304": "left out: widened from 36 neighbours, superseded by 46",
    "24_spotsweeper_8um_k48": "SpotSweeper size outcomes",
    "25_spotsweeper_16um_k48": "SpotSweeper size outcomes",
    "26_spotsweeper_32um_k48": "SpotSweeper size outcomes",
    "27_spotsweeper_64um_k48": "SpotSweeper size outcomes",
    "28_spotsweeper_128um_k48": "SpotSweeper size outcome for the largest class (512 um LQRs)",
    "40_spotsweeper_128um": "left out: 36-neighbour run of 28, superseded",
    "41_mad_128um": "left out: run on the 512 um subset only",
    "30_miqc_8um": "Discussion, avoid: miQC on 8 um bins at low depth",
    "31_miqc_64um": "left out: misses signal loss, no better than MAD on 64 um bins otherwise",
    "32_fixedumi_8um": "Discussion, avoid: a fixed UMI cutoff",
    "33_fixedumi_64um": "left out: a fixed UMI cutoff, covered by the Discussion's avoid step",
    "34_mad_64um_percovered": "left out: boundary correction that changed boundary loss little (129)",
    "35_mad_64um_excludepartial": "Discussion, check step: boundary correction and what it costs (129)",
    "44_spotsweeper_8um_k192": "SpotSweeper size outcomes, widened",
    "45_spotsweeper_8um_k768": "SpotSweeper size outcomes, widened",
    "46_spotsweeper_8um_k3072": "SpotSweeper size outcomes, widened",
}
COVERAGE_CONDITIONS = {
    "clean": "cost on undamaged tissue: the MAD bin and depth rules",
    **{f"{k}|{dm}": "SpotSweeper size outcomes" for k in KINDS for dm in LQRS},
    **{f"{k}|{VERY_LARGE}": "SpotSweeper size outcome for the largest class" for k in KINDS},
    **{f"{k}edge|{dm}": "edge decision and outcome" for k in KINDS for dm in LQRS},
}

def audit(S):
    cond = S.kind.where(S.kind == "clean", S.kind + "|" + S.diameter_um.fillna(0).astype(int).astype(str))
    arms, conds = sorted(S.arm.unique()), sorted(cond.unique())
    gaps = [a for a in arms if a not in COVERAGE_ARMS] + [c for c in conds if c not in COVERAGE_CONDITIONS]
    if gaps:
        raise SystemExit("not covered by the guideline or its exclusions: " + ", ".join(gaps))
    return {"arms": {a: COVERAGE_ARMS[a] for a in arms},
            "conditions": {c: COVERAGE_CONDITIONS[c] for c in conds}}

def med(s):
    s = pd.Series(s).dropna()
    return float(np.median(s)) if len(s) else None

def load():
    rows = []
    for sec, dep, cell in P.cells(GRID):
        f = P.scores(cell)
        f["section"], f["depth"], f["chem"] = sec, dep, P.chemistry(sec)
        rows.append(f)
    S = pd.concat(rows, ignore_index=True)
    S["kind"] = S["kind"].astype(str)
    return S

def numbers(S):
    G = {"generated": f"{pd.Timestamp.now(tz='Asia/Taipei'):%Y-%m-%d %H:%M}",
         "runs": int(S.groupby(["section", "depth"]).ngroups)}
    t1 = pd.read_csv(os.path.join(OUT, "table1_panel.csv"))
    G["panel"] = dict(sections=len(t1),
                      umi=[float(t1.native_mean_umi.min()), float(t1.native_mean_umi.max())],
                      mito=[float(t1.mito_baseline_pct.min()), float(t1.mito_baseline_pct.max())],
                      mito_n=int(t1.mito_baseline_pct.notna().sum()))

    def recall(arm, kind, dm, dep, ch):
        s = S[(S.arm == arm) & (S.kind == kind) & (S.diameter_um == dm) & (S.depth == dep) & (S.chem == ch)]
        return s.set_index("section").recall.dropna()

    def cost(arm, dep, ch):
        s = S[(S.arm == arm) & (S.condition == "clean") & (S.depth == dep) & (S.chem == ch)]
        return s.set_index("section").discard_frac.dropna()

    nb = {}
    for ch in P.CHEMS:
        for k in KINDS:
            for dm in LQRS:
                for dep in P.DEPTHS:
                    cur = {b: recall(SS[b], k, dm, dep, ch) for b in BINS}
                    common = sorted(set.intersection(*(set(c.index) for c in cur.values())))
                    if not common:
                        continue
                    m = {b: med(cur[b][common]) for b in BINS}
                    opt = max(BINS, key=lambda b: m[b])
                    e = dict(n=len(common), median=m, opt=opt)
                    w = {b: recall(WIDE[b], k, dm, dep, ch) for b in WIDE}
                    wc = sorted(set(common).intersection(*(set(c.index) for c in w.values())))
                    if wc:
                        e["widened"] = dict(n=len(wc), median={b: med(w[b][wc]) for b in WIDE},
                                            ad={b: med(cur[b][wc]) for b in WIDE})
                    nb[f"{ch}|{k}|{dm}|{dep}"] = e
    G["neighbourhood"] = nb

    nc = {}
    for ch in P.CHEMS:
        for dm in LQRS:
            for dep in P.DEPTHS:
                t = {}
                for b in BINS:
                    s = S[(S.arm == SS[b]) & (S.diameter_um == dm) & (S.depth == dep) & (S.chem == ch)
                          & S.kind.isin(KINDS)]
                    t[b] = s.pivot_table(index="section", columns="kind", values="recall")
                if not all(set(KINDS) <= set(x.columns) for x in t.values()):
                    continue
                common = sorted(set.intersection(*(set(x[list(KINDS)].dropna().index) for x in t.values())))
                if not common:
                    continue
                m = {b: med(t[b].loc[common, list(KINDS)].mean(axis=1)) for b in BINS}
                nc[f"{ch}|{dm}|{dep}"] = dict(n=len(common), median=m, opt=max(BINS, key=lambda b: m[b]))
    G["neighbourhood_combined"] = nc

    cand = {**{f"b{b}": SS[b] for b in BINS}, **{f"k{48 * (b // 8) ** 2}": WIDE[b] for b in WIDE}}
    wc = {}
    for dm in LQRS:
        for dep in ENDS:
            e = {}
            for k in KINDS:
                vals = {nm: recall(a, k, dm, dep, "probe") for nm, a in cand.items()}
                common = sorted(set.intersection(*(set(v.index) for v in vals.values())))
                e[k] = dict(n=len(common), median={nm: med(v[common]) for nm, v in vals.items()})
            wc[f"{dm}|{dep}"] = e
    G["widened_check"] = wc

    vl = {}
    arms_vl = {**SS, 128: "28_spotsweeper_128um_k48"}
    for ch in P.CHEMS:
        for dep in ENDS:
            per = {}
            for k in KINDS:
                vals = {b: recall(arms_vl[b], k, VERY_LARGE, dep, ch) for b in BINS_VL}
                common = sorted(set.intersection(*(set(v.index) for v in vals.values())))
                if not common:
                    continue
                per[k] = pd.DataFrame({b: vals[b][common] for b in BINS_VL})
                key = f"{k}|{dep}" if ch == "probe" else f"{ch}|{k}|{dep}"
                vl[key] = dict(n=len(common), median={b: med(per[k][b]) for b in BINS_VL})
            if len(per) == len(KINDS):
                both = sorted(set.intersection(*(set(x.index) for x in per.values())))
                m = {b: med(np.mean([per[k].loc[both, b] for k in KINDS], axis=0)) for b in BINS_VL}
                vl[f"combined|{ch}|{dep}"] = dict(n=len(both), median=m, opt=max(BINS_VL, key=lambda b: m[b]))
    G["very_large"] = vl

    bc = {}
    p = os.path.join(OUT, "boundary_corrections.csv")
    if os.path.exists(p):
        B = pd.read_csv(p)
        B = B[B.chem == "probe"]
        for dep in ENDS:
            x = B[B.depth == dep].assign(band=lambda d: d.edge_disc / d.edge_n)
            piv = x.pivot_table(index="section", columns="arm", values="band")
            core = [a for a in ("09_mad_32um", MAD[64], "34_mad_64um_percovered",
                                "35_mad_64um_excludepartial") if a in piv]
            pc = piv[core].dropna()
            for a in core:
                bc[f"band|{a}|{dep}"] = dict(n=len(pc), median=med(pc[a]))
            coarse = [a for a in (MAD[64], "41_mad_128um", SS[64], "28_spotsweeper_128um_k48") if a in piv]
            if any(a in coarse for a in ("41_mad_128um", "28_spotsweeper_128um_k48")):
                pk = piv[coarse].dropna()
                for a in coarse:
                    bc[f"band_coarse|{a}|{dep}"] = dict(n=len(pk), median=med(pk[a]))
    for a in (MAD[64], "35_mad_64um_excludepartial"):
        for dep in ENDS:
            r = recall(a, "dropout01edge", LQRS[-1], dep, "probe")
            bc[f"edge_recall|{a}|{dep}"] = dict(n=len(r), median=med(r))
    G["boundary"] = bc

    ed = {}
    for a in ("19_spatialartifacts_8um", "20_spatialartifacts_16um", *SS.values()):
        for dm in LQRS:
            for dep in ENDS:
                for k in ("dropout01edge", "mitoedge"):
                    r = recall(a, k, dm, dep, "probe")
                    ed[f"{a}|{k}|{dm}|{dep}"] = dict(n=len(r), median=med(r))
    G["edge"] = ed

    sw = {}
    for ch in P.CHEMS:
        for b in BINS:
            for dep in P.DEPTHS:
                c = cost(MAD[b], dep, ch)
                sw[f"{ch}|cost|{b}|{dep}"] = dict(n=len(c), median=med(c))
            a, z = cost(MAD[b], 100, ch), cost(MAD[b], 12, ch)
            j = a.index.intersection(z.index)
            rat = (z[j] / a[j].replace(0, np.nan)).dropna()
            sw[f"{ch}|ratio|{b}"] = dict(n=len(rat), median=med(rat))
            for k in KINDS:
                for dm in LQRS:
                    for dep in ENDS:
                        r = recall(MAD[b], k, dm, dep, ch)
                        sw[f"{ch}|recall|{b}|{k}|{dm}|{dep}"] = dict(n=len(r), median=med(r))
        for dep in P.DEPTHS:
            c = cost(SS[8], dep, ch)
            sw[f"{ch}|sscost|8|{dep}"] = dict(n=len(c), median=med(c))
    for b in BINS:
        k = n = 0
        for ch in P.CHEMS:
            a, z = cost(MAD[b], 100, ch), cost(MAD[b], 12, ch)
            j = a.index.intersection(z.index)
            n += len(j); k += int((z[j] > TOL * a[j]).sum())
        sw[f"all|rise|{b}"] = dict(n=n, k=k)
    G["section_wide"] = sw

    G["miqc_failed"] = {}
    for arm in ("30_miqc_8um", "31_miqc_64um"):
        fails = P.method_failures(GRID, arm, P.big512_runs(WK))
        if any(f for f, _ in fails.values()):
            G["miqc_failed"][arm] = {str(d): (f / n if n else None) for d, (f, n) in fails.items()}
    return G

def rules(G):
    sw = G["section_wide"]
    c = lambda ch, b, d: sw[f"{ch}|cost|{b}|{d}"]["median"]
    safe = next(b for b in BINS[1:] if all(
        max(c(ch, b, d) for d in P.DEPTHS) <= TOL * c(ch, b, 100) for ch in P.CHEMS))
    floor = 100
    for d in P.DEPTHS:
        if all(c(ch, 8, d) <= TOL * c(ch, safe, d) for ch in P.CHEMS):
            floor = d
        else:
            break
    return safe, floor


def settings(G):
    nc = G["neighbourhood_combined"]
    top, lower = P.DEPTHS[0], [d for d in P.DEPTHS[1:]]

    def tiers(ch, dm):
        full = nc[f"{ch}|{dm}|{top}"]["opt"]
        lo = [nc[f"{ch}|{dm}|{d}"]["opt"] for d in lower if f"{ch}|{dm}|{d}" in nc]
        most = max(lo.count(b) for b in lo)
        return full, min(b for b in lo if lo.count(b) == most)

    nbh = {dm: tiers("probe", dm) for dm in LQRS}
    agree = [dm for dm in LQRS if tiers(P.POLYA, dm) == nbh[dm]]
    cuts = [int(5 * round(math.sqrt(a * b) / 5)) for a, b in zip(LQRS, LQRS[1:])]
    safe, floor = rules(G)

    wc = G["widened_check"]

    def widen(dm, dep, base):
        e = wc.get(f"{dm}|{dep}", {})
        kinds = [k for k in KINDS if e.get(k, {}).get("n", 0) >= 1]
        prov = any(e.get(k, {}).get("n", 0) < WIDE_MIN for k in KINDS)
        if not kinds:
            return None, prov
        score = lambda nm: float(np.mean([e[k]["median"][nm] for k in kinds]))
        ks = [nm for nm in e[kinds[0]]["median"] if nm.startswith("k")]
        best = max(score(nm) for nm in ks)
        if best <= score(f"b{base}") + 0.01:
            return None, prov
        return min((int(nm[1:]) for nm in ks if score(nm) >= best - 0.01)), prov

    rec = {}
    for dm in LQRS:
        full, low = nbh[dm]
        for tier, dep, base in (("hi", top, full), ("lo", ENDS[-1], low)):
            k, prov = widen(dm, dep, base)
            rec[(dm, tier)] = dict(bin=8, k=k, prov=prov) if k else dict(bin=base, k=48, prov=False)

    ed = G["edge"]
    sa = {}
    for tier, dep in (("hi", top), ("lo", ENDS[-1])):
        sc = {b: np.mean([ed[f"{a}|dropout01edge|{dm}|{dep}"]["median"] or 0 for dm in LQRS[1:]])
              for b, a in ((8, "19_spatialartifacts_8um"), (16, "20_spatialartifacts_16um"))}
        sa[tier] = max(sc, key=sc.get)

    vl = G.get("very_large", {})
    sizes, extra = list(LQRS), []
    if all(f"combined|probe|{d}" in vl for d in ENDS):
        sizes.append(VERY_LARGE)
        hi_b, lo_b = (vl[f"combined|probe|{d}"]["opt"] for d in ENDS)
        extra.append(dict(hi=dict(bin=hi_b, k=48, prov=False), lo=dict(bin=lo_b, k=48, prov=False), note=None))
    cuts = [int(5 * round(math.sqrt(a * b) / 5)) for a, b in zip(sizes, sizes[1:])]
    labels = ([f"< {cuts[0]} µm"] + [f"{a}–{b} µm" for a, b in zip(cuts, cuts[1:])]
              + [f"> {cuts[-1]} µm"])
    classes = []
    for dm, label in zip(LQRS, labels):
        hi, lo = rec[(dm, "hi")], rec[(dm, "lo")]
        note = "more than 48 neighbours: provisional, see Fig. S9" if hi.get("prov") or lo.get("prov") else None
        classes.append(dict(label=label, hi=hi, lo=lo, note=note))
    for e, label in zip(extra, labels[len(LQRS):]):
        classes.append(dict(label=label, **e))
    return dict(nbh=nbh, agree=agree, cuts=cuts, safe=safe, floor=floor, top=top, rec=rec, sa=sa,
                classes=classes, sizes=sizes)

def lib():
    spec = importlib.util.spec_from_file_location(
        "methods_diagrams", os.path.join(WK, "scripts", "117_methods_diagrams.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def draw(G):
    D = lib()
    D.W = 820
    C, text, rect, path = D.C, D.text, D.rect, D.path
    st = settings(G)
    top, cuts, rec, sa = st["top"], st["cuts"], st["rec"], st["sa"]
    safe, floor = st["safe"], st["floor"]
    b = []

    def card(x, y, w, h, title, node=""):
        return [rect(x, y, w, h, C["paper"], rx=6),
                rect(x, y, w, h, C["paper"], C["ink"], 1, rx=6, extra=f' data-node="{node}"'),
                text(x + 16, y + h / 2 + 5, title, 15, weight=700, node=node)]

    def outcome(x, y, w, h, tool, rows, notes=(), node=""):
        out = [rect(x, y, w, h, C["paper"], rx=6),
               rect(x, y, w, h, C["paper"], C["ink"], 1.4, rx=6, extra=f' data-node="{node}"'),
               text(x + 12, y + 24, tool, 15, weight=700, node=node)]
        yy = y + 24
        for lab, val in rows:
            yy += 19
            out += [text(x + 12, yy, lab, 11, C["muted"], node=node),
                    text(x + 108, yy, val, 12, weight=700, node=node)]
        for ln in notes:
            yy += 16
            out.append(text(x + 12, yy, ln, 11, C["muted"], italic=True, node=node))
        return out

    def diamond(cx, cy, hw, hh, lines, node):
        pts = f"{cx},{cy - hh} {cx + hw},{cy} {cx},{cy + hh} {cx - hw},{cy}"
        out = [f'<polygon points="{pts}" fill="{C["paper"]}" stroke="{C["ink"]}" stroke-width="1.2" '
               f'data-node="{node}"/>']
        y0 = cy - 8 * (len(lines) - 1) + 5
        out += [text(cx, y0 + 16 * i, s, 13, weight=700, anchor="middle", node=node)
                for i, s in enumerate(lines)]
        return out

    def lab(x, y, s, anchor="start"):
        w = 6.8 * len(s) + 8
        x0 = x - 4 if anchor == "start" else x - w / 2
        return [rect(x0, y - 12, w, 16, C["paper"]),
                text(x, y, s, 11, C["muted"], anchor=anchor, fam=D.MONO)]

    def setting(r):
        return f"{r['bin']} µm bins" + ("" if r["bin"] == 8 else " + A&D") + f", {r['k']} neighbours"

    hi, lo = f"≥ {top} UMIs/bin", f"< {top} UMIs/bin"
    TX, TW = 440, 340
    RX, RW = 544, 236
    cls = st["classes"]
    same = lambda c: (c["hi"]["bin"], c["hi"]["k"]) == (c["lo"]["bin"], c["lo"]["k"])
    hs = [(48 if same(c) else 68) + (24 if c["note"] else 0) for c in cls]
    S1, S2, CH = 20, 84, 44
    D1 = (230, S2 + CH + 20 + 46)
    DE = (230, D1[1] + 46 + 82 + 44)
    T_MAD = (RX, D1[1] - 26, RW, 76); T_EDGE = (RX, DE[1] - 26, RW, 76)
    rejoin = T_EDGE[1] + T_EDGE[3] + 16
    ys, y = [], rejoin + 70
    for h in hs:
        ys.append(y); y += h + 20
    cy = [y_ + 26 for y_ in ys]
    DS = (230, round((cy[0] + cy[-1]) / 2))
    y_chk = ys[-1] + hs[-1] + 32

    b += [path([(230, S1 + CH), (230, S2)]), path([(230, S2 + CH), (230, D1[1] - 46)]),
          path([(360, D1[1]), (RX, D1[1])]),
          path([(230, D1[1] + 46), (230, DE[1] - 44)]),
          path([(350, DE[1]), (RX, DE[1])]),
          path([(230, DE[1] + 44), (230, DS[1] - 48)]),
          path([(RX + RW / 2, T_EDGE[1] + T_EDGE[3]), (RX + RW / 2, rejoin), (230, rejoin)], marker=None),
          path([(340, DS[1]), (356, DS[1])], marker=None),
          path([(356, cy[0]), (356, cy[-1])], marker=None)]
    b += [path([(356, c), (TX, c)]) for c in cy]
    ends = [D1[1]] + cy
    b += [path([(780, e), (800, e)], marker=None) for e in ends]
    b += [path([(800, ends[0]), (800, y_chk + CH / 2), (780, y_chk + CH / 2)])]
    b += [f'<circle cx="{x}" cy="{y}" r="3" fill="{C["muted"]}"/>'
          for x, y in [(230, rejoin), (356, DS[1])] + [(356, c) for c in cy[1:-1]] + [(800, e) for e in ends[1:]]]
    b += lab(452, D1[1] - 8, "NO / NOT SURE", "middle") + lab(238, (D1[1] + DE[1]) / 2, "YES")
    b += lab(447, DE[1] - 8, "YES", "middle") + lab(238, rejoin + 18, "NO")
    for c, cl in zip(cy, cls):
        b += lab(398, c - 8, cl["label"], "middle")

    b += card(20, S1, 760, CH, "Compute QC metrics on the 8 µm bins", node="s1")
    b += card(20, S2, 760, CH, "Inspect the tissue image", node="s2")

    b += diamond(*D1, 130, 46, ("Damage visible and", "its size known?"), "d-seen")
    b += diamond(*DE, 120, 44, ("At the tissue", "edge?"), "d-edge")
    b += diamond(*DS, 110, 48, ("How large are", "the patches?"), "d-size")

    b += outcome(*T_MAD, "MAD threshold, scuttle", [(hi, "8 µm bins"), (lo, f"{safe} µm bins + A&D")],
                 node="t-mad")
    b += outcome(*T_EDGE, "Also SpatialArtifacts", [(hi, f"{sa['hi']} µm bins" + ("" if sa["hi"] == 8 else " + A&D")),
                                                    (lo, f"{sa['lo']} µm bins" + ("" if sa["lo"] == 8 else " + A&D"))],
                 node="t-edge")
    for i, cl in enumerate(cls):
        rows = ([("any depth", setting(cl["hi"]))] if same(cl)
                else [(hi, setting(cl["hi"])), (lo, setting(cl["lo"]))])
        b += outcome(TX, ys[i], TW, hs[i], "SpotSweeper", rows, (cl["note"],) if cl["note"] else (),
                     f"t-size{i}")

    b += card(20, y_chk, 760, CH, "Check what was removed", node="s3")
    H = y_chk + CH + 20

    wc = G["widened_check"][f"{LQRS[0]}|{top}"]
    cap = (f"<b>How to run QC on a Visium HD section.</b> Each outcome names the package, the bin "
           f"size and the number of neighbours to run it with. Settings are chosen from the benchmark "
           f"by the code that draws this figure: for SpotSweeper, the bin at which recall averaged over signal "
           f"loss and mitochondrial elevation peaked (probe-based sections), replaced by a wider "
           f"neighbourhood on 8 µm bins where that did better on the sections where both were run "
           f"(signal loss n = {wc['dropout01']['n']}, mitochondrial elevation n = {wc['mito']['n']}); "
           + (f"the largest band from the sections that hold 512 µm LQRs (n = "
              f"{G['very_large'][f'combined|probe|{top}']['n']} at {top} UMIs per bin, "
              f"{G['very_large'][f'combined|probe|{ENDS[-1]}']['n']} at {ENDS[-1]}); "
              if VERY_LARGE in st["sizes"] else "") +
           f"size bands at the geometric midpoints between the diameters tested, "
           f"{', '.join(map(str, st['sizes'][:-1]))} and {st['sizes'][-1]} µm; for SpatialArtifacts, the bin that "
           f"found more low-count damage at the tissue edge; for the MAD threshold, the finest bin "
           f"whose cost on undamaged tissue did not grow as depth fell. Drawn from {G['runs']} scored "
           f"runs, {G['generated']}.")
    desc = ("Flowchart: compute QC metrics on 8 µm bins and inspect the image. If no damage is seen, "
            "run a MAD threshold; if damage is seen, add SpatialArtifacts at the tissue edge and run "
            "SpotSweeper with a bin size and number of neighbours set by the size of the damage and the depth; "
            "then check what was removed.")
    return D, D.page("guideline", "Discussion · guideline", "How to run QC on a Visium HD section",
                     desc, H, "\n".join(b), cap)

if __name__ == "__main__":
    S = load()
    cover = audit(S)
    G = json.loads(json.dumps(numbers(S)))
    G["coverage"] = cover
    print(f"coverage: {len(cover['arms'])} arms, {len(cover['conditions'])} conditions, each used or excused")
    p = os.path.join(OUT, "guideline_numbers.json")
    json.dump(G, open(p, "w"), indent=1)
    print("wrote", p)
    st = settings(G)
    bounds = [0] + st["cuts"] + [None]
    tool = {"generated": G["generated"], "depth_cut_umi_per_bin": st["top"],
            "classes": [dict(label=c["label"], min_um=bounds[i], max_um=bounds[i + 1],
                             hi=dict(c["hi"]), lo=dict(c["lo"]))
                        for i, c in enumerate(st["classes"])],
            "mad": {"hi": 8, "lo": st["safe"]}, "spatialartifacts": st["sa"]}
    p = os.path.join(OUT, "guideline_settings.json")
    json.dump(tool, open(p, "w"), indent=1)
    print("wrote", p)
    D, html = draw(G)
    base = os.path.join(OUT, "fig_guideline")
    open(base + ".html", "w").write(html)
    i, j = html.index("<svg "), html.index("</svg>") + len("</svg>")
    open(base + ".svg", "w").write('<?xml version="1.0" encoding="UTF-8"?>\n' + html[i:j].replace(
        'role="img"', f'width="{D.W}" role="img"', 1) + "\n")
    import shutil, subprocess
    if shutil.which("rsvg-convert"):
        subprocess.run(["rsvg-convert", "-z", "3", "-o", base + ".png", base + ".svg"], check=True)
        subprocess.run(["rsvg-convert", "-f", "pdf", "-o", base + ".pdf", base + ".svg"], check=True)
    print("wrote", base + ".{html,svg,png,pdf}")

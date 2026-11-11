#!/usr/bin/env python3
import argparse, glob, importlib.util, json, os, re, sys
import numpy as np, pandas as pd

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GRID = os.path.join(WK, "benchmark_grid")
FIG = os.path.join(GRID, "figures")
SETTINGS = os.path.join(FIG, "guideline_settings.json")

def arm(results, pattern):
    hits = [a for a in os.listdir(results) if re.fullmatch(pattern, a)]
    if len(hits) > 1:
        raise SystemExit(f"{pattern} matches {hits}")
    return hits[0] if hits else pattern

def pairs(results, deep, shallow):
    cfg = json.load(open(SETTINGS))
    assert deep >= cfg["depth_cut_umi_per_bin"] > shallow
    ss = lambda b, k: arm(results, rf"\d+_spotsweeper_{b}um_k{k}")
    mad = lambda b: arm(results, rf"\d+_mad_{b}um")
    sa = lambda b: arm(results, rf"\d+_spatialartifacts_{b}um")
    out = [("no_filter", [], []),
           ("mad_64um_both", [mad(64)], [mad(64)]),
           ("procedure_section_wide", [mad(cfg["mad"]["hi"])], [mad(cfg["mad"]["lo"])])]
    s = cfg["spatialartifacts"]
    for c in cfg["classes"]:
        hi, lo = [ss(c["hi"]["bin"], c["hi"]["k"])], [ss(c["lo"]["bin"], c["lo"]["k"])]
        if c["min_um"] >= s["min_um"]:
            hi.append(sa(s["hi"]))
            lo.append(sa(s["lo"]))
        name = f"procedure_{c['min_um']}_{c['max_um'] if c['max_um'] else 'up'}um"
        out.append((name, hi, lo))
    return out

def removed(cell, arms):
    drop = set()
    for a in arms:
        f = os.path.join(cell, "results", a, "clean.txt")
        if not os.path.exists(f):
            return None
        drop |= set(open(f).read().split())
    return drop

def run(a):
    spec = importlib.util.spec_from_file_location("ds94", os.path.join(WK, "scripts", "94_downstream_binned.py"))
    ds = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ds)
    from sklearn.metrics import adjusted_rand_score

    cells = {d: os.path.join(GRID, f"{a.section}_d{d}") for d in (a.deep, a.shallow)}
    raw = {d: ds.read_bins(c) for d, c in cells.items()}
    for d, x in raw.items():
        print(f"{a.section} d{d}: {x.n_obs} bins of 8 µm, {x.n_vars} genes", flush=True)

    rows = []
    for name, hi, lo in pairs(os.path.join(cells[a.deep], "results"), a.deep, a.shallow):
        arms = {a.deep: hi, a.shallow: lo}
        drop = {d: removed(cells[d], arms[d]) for d in cells}
        if any(v is None for v in drop.values()):
            print(f"  {name}: result missing for {arms}, skipped", flush=True)
            continue
        lab, sub = {}, {}
        for dep in cells:
            keep = [b for b in raw[dep].obs_names if b not in drop[dep]]
            sub[dep] = ds.pool(raw[dep][keep].copy(), a.bin)
            if sub[dep].n_obs < 500:
                print(f"  {name} d{dep}: only {sub[dep].n_obs} bins, skipped", flush=True)
                break
            lab[dep] = ds.cluster(sub[dep], seed=0)
            print(f"  {name} d{dep}: {'+'.join(arms[dep]) or 'none'}, {len(keep)} of {raw[dep].n_obs} "
                  f"8 µm bins kept, {sub[dep].n_obs} bins of {a.bin} µm, {lab[dep].nunique()} clusters", flush=True)
        if len(lab) != 2:
            continue
        base = dict(section=a.section, arm=name, bin_um=a.bin,
                    arms_a="+".join(hi) or "none", arms_b="+".join(lo) or "none")
        reseed = ds.cluster(sub[a.deep], seed=1)
        sh = lab[a.deep].index.intersection(reseed.index)
        rows.append(dict(base, comparison="reseed", depth_a=a.deep, depth_b=a.deep, n_shared=len(sh),
                         ari=adjusted_rand_score(lab[a.deep][sh], reseed[sh]),
                         n_clusters_a=lab[a.deep].nunique(), n_clusters_b=reseed.nunique(),
                         kept_a=sub[a.deep].n_obs, kept_b=sub[a.deep].n_obs))
        sh = lab[a.deep].index.intersection(lab[a.shallow].index)
        rows.append(dict(base, comparison="across_depth", depth_a=a.deep, depth_b=a.shallow, n_shared=len(sh),
                         ari=adjusted_rand_score(lab[a.deep][sh], lab[a.shallow][sh]),
                         n_clusters_a=lab[a.deep].nunique(), n_clusters_b=lab[a.shallow].nunique(),
                         kept_a=sub[a.deep].n_obs, kept_b=sub[a.shallow].n_obs))
        print(f"  -> {name} ARI across depth {rows[-1]['ari']:.3f} on {len(sh)} shared bins, "
              f"reseed {rows[-2]['ari']:.3f}", flush=True)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    tmp = a.out + ".part"
    pd.DataFrame(rows).to_csv(tmp, index=False)
    os.replace(tmp, a.out)
    print("wrote", len(rows), "rows to", a.out)

def sign_p(k, n):
    from scipy.stats import binomtest
    return binomtest(k, n, 0.5).pvalue if n else float("nan")

def report():
    sys.path.insert(0, os.path.join(WK, "scripts"))
    import panel as P
    files = sorted(glob.glob(os.path.join(FIG, "downstream_guideline", "*.csv")))
    D = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    D.to_csv(os.path.join(FIG, "downstream_guideline.csv"), index=False)
    D["chem"] = D.section.map(P.chemistry)
    X = D[D.comparison == "across_depth"].pivot(index="section", columns="arm", values="ari")
    R = D[D.comparison == "reseed"].pivot(index="section", columns="arm", values="ari")
    chem = D.drop_duplicates("section").set_index("section").chem
    print(f"{len(files)} sections")

    old = pd.read_csv(os.path.join(FIG, "downstream_binned.csv"))
    old = old[old.comparison == "across_depth"].pivot(index="section", columns="arm", values="ari")
    for new, ref in (("no_filter", "01_no_qc"), ("mad_64um_both", "10_mad_64um")):
        d = (X[new] - old.reindex(X.index)[ref]).abs()
        print(f"check {new} against {ref} in downstream_binned.csv: largest difference {d.max():.2e} over {d.notna().sum()} sections")

    from scipy.stats import wilcoxon
    for c in ("probe", "3prime"):
        idx = chem[chem == c].index
        print(f"\n{c} sections (n = {len(idx)})")
        print(f"{'arm':30} {'n':>3} {'ARI':>6} {'reseed':>6} {'diff vs none':>12} {'lower':>7} {'sign p':>8} {'wilcoxon p':>10}")
        for a in X.columns:
            x = X.loc[idx, a].dropna()
            row = f"{a:30} {len(x):3d} {x.median():6.3f} {R.loc[idx, a].median():6.3f}"
            if a != "no_filter":
                d = (X.loc[x.index, a] - X.loc[x.index, "no_filter"]).dropna()
                k, n = int((d < 0).sum()), int((d != 0).sum())
                w = wilcoxon(d).pvalue if n else float("nan")
                row += f" {d.median():12.3f} {k:3d}/{n:<3d} {sign_p(k, n):8.4f} {w:10.4f}"
            print(row)
        for a in [c for c in X.columns if c.startswith("procedure_")]:
            d = (X.loc[idx, a] - X.loc[idx, "mad_64um_both"]).dropna()
            k, n = int((d < 0).sum()), int((d != 0).sum())
            print(f"  {a} against mad_64um_both: median {d.median():+.3f}, lower in {k}/{n}, sign p {sign_p(k, n):.4f}")

if __name__ == "__main__":
    if sys.argv[1:2] == ["report"]:
        report()
    else:
        ap = argparse.ArgumentParser()
        ap.add_argument("--section", required=True)
        ap.add_argument("--deep", type=int, default=100)
        ap.add_argument("--shallow", type=int, default=12)
        ap.add_argument("--bin", type=int, default=32)
        ap.add_argument("--out", required=True)
        run(ap.parse_args())

import glob, json, os, shutil, subprocess, sys
import numpy as np, pandas as pd
sys.path.insert(0, "/work/scripts")
import importlib.util
spec = importlib.util.spec_from_file_location("inj", "/work/scripts/10_inject_artifacts.py")
inj = importlib.util.module_from_spec(spec); spec.loader.exec_module(inj)

run, extra = sys.argv[1], sys.argv[2:]
R = f"/work/benchmark_grid/{run}"
TMP = f"/work/deposit/rebuild_tmp/{run}"
shutil.rmtree(TMP, ignore_errors=True)
on_disk = {os.path.basename(os.path.dirname(g)) for g in glob.glob(f"{R}/*/ground_truth_008um.parquet")}
written, bad = set(), 0
passed, differs = {}, {}
only = os.environ.get("ONLY_GIVEN") == "1"
for pf in (extra if only else sorted(glob.glob(f"{R}/injection_params*.json")) + extra):
    p = json.load(open(pf))
    out = os.path.join(TMP, os.path.basename(pf)[:-5])
    cmd = [sys.executable, "/work/scripts/10_inject_artifacts.py", "--src", "/src", "--out", out,
           "--diameters", *map(str, p["diameters"]), "--n-regions", str(p["n_regions"]),
           "--keep-p", str(p["keep_p"]), "--mito-target-frac", repr(p["mito_target_frac"]),
           "--seed", str(p["seed"]), "--target-umi-per-bin", repr(p["target_umi_per_bin"]),
           "--kinds", *p["kinds"], "--location", p.get("location") or "interior"]
    for k, flag in (("mito_odds_ratio", "--mito-odds-ratio"), ("mito_baseline", "--mito-baseline"),
                    ("kind_label", "--kind-label"), ("record_suffix", "--record-suffix")):
        if p.get(k) not in (None, ""):
            cmd += [flag, str(p[k])]
    if p.get("with_clean"):
        cmd.append("--with-clean")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode:
        print(f"FAIL {os.path.basename(pf)}: injector error\n{r.stderr[-600:]}"); bad = 1; continue
    for c in sorted(d for d in os.listdir(out) if os.path.isdir(os.path.join(out, d))):
        if c not in on_disk:
            print(f"SUPERSEDED {c} (from {os.path.basename(pf)}): no longer on disk", flush=True)
            continue
        written.add(c)
        a, b = os.path.join(out, c), os.path.join(R, c)
        m1, c1, _ = inj.read_mex(f"{a}/square_008um/filtered_feature_bc_matrix")
        m2, c2, _ = inj.read_mex(f"{b}/square_008um/filtered_feature_bc_matrix")
        g1 = pd.read_parquet(f"{a}/ground_truth_008um.parquet"); g2 = pd.read_parquet(f"{b}/ground_truth_008um.parquet")
        p1 = pd.read_parquet(f"{a}/square_008um/spatial/tissue_positions.parquet")
        p2 = pd.read_parquet(f"{b}/square_008um/spatial/tissue_positions.parquet")
        ok = {"barcodes": list(c1) == list(c2),
              "matrix": m1.shape == m2.shape and (m1 != m2).nnz == 0,
              "ground_truth": g1.equals(g2),
              "positions": p1.equals(p2)}
        (passed if all(ok.values()) else differs).setdefault(c, []).append(os.path.basename(pf))
        print(("PASS " if all(ok.values()) else "DIFFERS ") + f"{c} (from {os.path.basename(pf)}): "
              + ", ".join(f"{k} {'same' if v else 'DIFFERENT'}" for k, v in ok.items()), flush=True)
    shutil.rmtree(out, ignore_errors=True)
for c in sorted(differs):
    if c in passed:
        print(f"SUPERSEDED {c}: {', '.join(differs[c])} replaced by {', '.join(passed[c])}")
checked = (passed.keys() | differs.keys()) if only else on_disk
fail = sorted(c for c in checked if c not in passed)
for c in fail:
    print(f"NOT REBUILT {c}: " + (f"differs from {', '.join(differs[c])}" if c in differs
                                  else "no recorded injector call wrote it"))
shutil.rmtree(TMP, ignore_errors=True)
print(f"ALL PASS ({len(passed)} conditions)" if not bad and not fail else f"NOT ALL REBUILT ({len(fail)} of {len(checked)})")

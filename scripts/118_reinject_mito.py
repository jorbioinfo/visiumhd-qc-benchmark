import argparse, json, os, subprocess, sys
import pandas as pd

FAMILIES = {
    "core": (["--diameters", "32", "64", "256"], ["mito_32um", "mito_64um", "mito_256um"]),
    "big":  (["--diameters", "512", "--n-regions", "8"], ["mito_512um"]),
    "edge": (["--diameters", "32", "64", "256", "--location", "edge"],
             ["mitoedge_32um", "mitoedge_64um", "mitoedge_256um"]),
}

ap = argparse.ArgumentParser()
ap.add_argument("--tag", required=True)
ap.add_argument("--depth", required=True)
ap.add_argument("--families", required=True)
ap.add_argument("--archive", required=True)
a = ap.parse_args()
out = f"/work/benchmark_grid/{a.tag}"
report = {"tag": a.tag, "families": {}, "ok": True}
for fam in a.families.split(","):
    extra, conds = FAMILIES[fam]
    cmd = ["/opt/pyenv/bin/python", "/work/scripts/10_inject_artifacts.py", "--src", "/src",
           "--out", out, "--target-umi-per-bin", a.depth, "--kinds", "mito",
           "--mito-odds-ratio", "10", "--mito-baseline", "native-mean",
           "--record-suffix", f"_mito_relative_{fam}"] + extra
    r = subprocess.run(cmd, capture_output=True, text=True)
    log = r.stdout + r.stderr
    target = [l for l in log.splitlines() if "odds ratio" in l]
    fr = {"exit": r.returncode, "target": target[-1].strip() if target else None, "conds": {}}
    if r.returncode != 0:
        fr["error"] = log[-600:]; report["ok"] = False
    for c in conds:
        new_f = f"{out}/{c}/ground_truth_008um.parquet"
        old_f = f"{a.archive}/{c}/ground_truth_008um.parquet"
        if not (os.path.exists(new_f) and os.path.exists(old_f)):
            fr["conds"][c] = "missing"; report["ok"] = False; continue
        n, o = pd.read_parquet(new_f), pd.read_parquet(old_f)
        same = set(n.barcode[n.is_artifact]) == set(o.barcode[o.is_artifact])
        fr["conds"][c] = f"regions identical ({int(n.is_artifact.sum())} bins)" if same else \
                         f"REGIONS DIFFER new {int(n.is_artifact.sum())} old {int(o.is_artifact.sum())}"
        report["ok"] &= same
    report["families"][fam] = fr
    with open(f"{out}/reinject_mito_{fam}.log", "w") as fh:
        fh.write(log)
print(json.dumps(report))
sys.exit(0 if report["ok"] else 3)

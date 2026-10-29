import argparse, glob, json, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARM_ORDER = [
    "01_no_qc", "02_mad_8um", "03_spotsweeper_8um", "04_spotsweeper_16um",
    "05_spotsweeper_32um", "06_spotsweeper_64um", "07_spatialartifacts",
]
ARM_LABEL = {
    "01_no_qc": "no QC", "02_mad_8um": "global MAD\n8 µm",
    "03_spotsweeper_8um": "SpotSweeper\n8 µm", "04_spotsweeper_16um": "SpotSweeper\n16 µm",
    "05_spotsweeper_32um": "SpotSweeper\n32 µm", "06_spotsweeper_64um": "SpotSweeper\n64 µm\n(current)",
    "07_spatialartifacts": "SpatialArtifacts",
}

def score(truth, pred):
    tp = len(truth & pred); fp = len(pred - truth); fn = len(truth - pred)
    prec = tp / (tp + fp) if tp + fp else np.nan
    rec = tp / (tp + fn) if tp + fn else np.nan
    f1 = 2 * prec * rec / (prec + rec) if tp and prec + rec else 0.0
    return dict(tp=tp, fp=fp, fn=fn, precision=prec, recall=rec, f1=f1)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--benchmark", required=True, help="the run directory")
    ap.add_argument("--results", required=True, help="the run's results directory")
    ap.add_argument("--out", required=True, help="where to write the scores")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows = []
    for gt_path in sorted(glob.glob(f"{a.benchmark}/*/ground_truth_008um.parquet")):
        cond = os.path.basename(os.path.dirname(gt_path))
        gt = pd.read_parquet(gt_path)
        truth = set(gt.barcode[gt.is_artifact])

        pos_path = os.path.join(os.path.dirname(gt_path), "square_008um",
                                "spatial", "tissue_positions.parquet")
        universe = set(gt.barcode)
        if os.path.exists(pos_path):
            pos = pd.read_parquet(pos_path)
            if "in_tissue" in pos.columns:
                universe = set(pos.barcode[pos.in_tissue == 1])
                truth &= universe
        kind = "clean" if cond == "clean" else cond.split("_")[0]
        diam = np.nan if cond == "clean" else int(cond.split("_")[1].replace("um", ""))

        for arm_dir in sorted(glob.glob(f"{a.results}/*")):
            arm = os.path.basename(arm_dir)
            f = os.path.join(arm_dir, f"{cond}.txt")
            if not os.path.exists(f):
                continue
            with open(f) as fh:
                pred = {l.strip() for l in fh if l.strip()}
            unknown = pred - universe
            if unknown:
                print(f"warning: {arm}/{cond} listed {len(unknown)} barcode(s) "
                      f"absent from the 8 µm table; ignoring them")
                pred &= universe
            r = score(truth, pred)
            r.update(arm=arm, condition=cond, kind=kind, diameter_um=diam,
                     n_pred=len(pred), n_truth=len(truth), n_total=len(universe),
                     discard_frac=len(pred) / len(universe))
            rows.append(r)

    if not rows:
        print(f"no arm results found under {a.results}/<arm>/<condition>.txt")
        return
    df = pd.DataFrame(rows)
    df.to_csv(f"{a.out}/scores.csv", index=False)
    print(df.to_string(index=False))

    arms = [x for x in ARM_ORDER if x in set(df.arm)] + \
           sorted(set(df.arm) - set(ARM_ORDER))
    cmap = plt.get_cmap("viridis")
    colors = {a_: cmap(i / max(len(arms) - 1, 1)) for i, a_ in enumerate(arms)}

    for kind in sorted(set(df.kind) - {"clean"}):
        d = df[df.kind == kind]
        diams = sorted(d.diameter_um.dropna().unique())
        fig, ax = plt.subplots(figsize=(1.6 * len(diams) + 3, 4.2))
        w = 0.8 / len(arms)
        for i, arm in enumerate(arms):
            vals = [d[(d.arm == arm) & (d.diameter_um == x)].f1.mean() for x in diams]
            ax.bar(np.arange(len(diams)) + i * w - 0.4 + w / 2, vals, w,
                   label=ARM_LABEL.get(arm, arm), color=colors[arm])
        ax.set_xticks(range(len(diams)))
        ax.set_xticklabels([f"{int(x)} µm" for x in diams])
        ax.set_xlabel("artifact diameter"); ax.set_ylabel("F1"); ax.set_ylim(0, 1)
        titles = {"dropout": "signal dropout, 20% of counts retained",
                  "dropout05": "signal dropout, 5% of counts retained",
                  "mito": "elevated mitochondrial fraction"}
        ax.set_title(f"Detection of injected artifacts — {titles.get(kind, kind)}")
        ax.legend(fontsize=7, ncol=2, frameon=False)
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout(); fig.savefig(f"{a.out}/detection_{kind}.png", dpi=160)
        print(f"wrote {a.out}/detection_{kind}.png")

    c = df[df.condition == "clean"]
    if len(c):
        fig, ax = plt.subplots(figsize=(1.3 * len(arms) + 2, 4.2))
        vals = [c[c.arm == arm].discard_frac.mean() * 100 for arm in arms]
        ax.bar(range(len(arms)), vals, color=[colors[a_] for a_ in arms])
        for i, v in enumerate(vals):
            if np.isfinite(v):
                ax.text(i, v, f"{v:.2f}%", ha="center", va="bottom", fontsize=8)
        ax.set_xticks(range(len(arms)))
        ax.set_xticklabels([ARM_LABEL.get(a_, a_) for a_ in arms], fontsize=8)
        ax.set_ylabel("good tissue discarded (%)")
        ax.set_title("False positives on uncorrupted data")
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout(); fig.savefig(f"{a.out}/false_positive.png", dpi=160)
        print(f"wrote {a.out}/false_positive.png")

if __name__ == "__main__":
    main()

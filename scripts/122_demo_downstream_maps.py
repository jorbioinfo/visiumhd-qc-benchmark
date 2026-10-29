import argparse, importlib.util, os, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score

_DS94 = []

def ds94():
    if not _DS94:
        spec = importlib.util.spec_from_file_location("ds94", os.path.join(os.path.dirname(__file__),
                                                                           "94_downstream_binned.py"))
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); _DS94.append(m)
    return _DS94[0]

ARMS = [("01_no_qc", "No QC"),
        ("02_mad_8um", "MAD, native 8 µm"),
        ("10_mad_64um", "A&D: MAD on 64 µm bins"),
        ("24_spotsweeper_8um_k48", "Neighbourhood test\n(SpotSweeper, k = 48)")]
BC = re.compile(r"s_(\d+)um_(\d+)_(\d+)")
ROW = {"01_no_qc": "No filter", "02_mad_8um": "MAD,\nnative 8 µm", "10_mad_64um": "MAD + A&D,\n64 µm bins",
       "24_spotsweeper_8um_k48": "SpotSweeper,\nnative 8 µm"}

def square_of(barcodes, bin_um=32):
    step = bin_um // 8
    return np.array([f"r{int(m.group(2)) // step}_c{int(m.group(3)) // step}"
                     for m in map(BC.match, barcodes)], dtype=object)

def rc_of(names):
    rc = np.array([[int(x) for x in re.findall(r"\d+", n)] for n in names])
    return rc[:, 0], rc[:, 1]

def match(lab, ref, tag=""):
    sh = lab.index.intersection(ref.index)
    ct = pd.crosstab(lab[sh], ref[sh])
    r, c = linear_sum_assignment(-ct.to_numpy())
    m = {ct.index[i]: ct.columns[j] for i, j in zip(r, c)}
    extra = [x for x in lab.unique() if x not in m]
    m.update({x: f"new{tag}_{k}" for k, x in enumerate(sorted(extra))})
    return lab.map(m)

def removed_fraction(names8, keep8, bin_um=32):
    sq = square_of(names8, bin_um)
    tot = pd.Series(1, index=sq).groupby(level=0).sum()
    kept = pd.Series(keep8.astype(int), index=sq).groupby(level=0).sum()
    return 1 - kept / tot

def removal(section, grid):
    rem, frac = {}, {}
    for d in (100, 12):
        cell = os.path.join(grid, f"{section}_d{d}")
        bc = pd.read_csv(os.path.join(cell, "clean", "square_008um", "filtered_feature_bc_matrix",
                                      "barcodes.tsv.gz"), header=None)[0].to_numpy()
        for arm, _ in ARMS:
            f = os.path.join(cell, "results", arm, "clean.txt")
            drop = set(open(f).read().split()) if os.path.exists(f) else set()
            keep8 = np.array([b not in drop for b in bc])
            rem[arm, d], frac[arm, d] = removed_fraction(bc, keep8), 1 - keep8.mean()
    return rem, frac

def replot(section, grid, out):
    df = pd.read_csv(os.path.join(out, f"{section}_labels.csv.gz"), dtype={"label": str})
    df["agree"] = df["agree"].map({True: True, False: False, "True": True, "False": False})
    sm = pd.read_csv(os.path.join(out, "summary.csv")).query("section == @section").set_index("arm")
    rem, frac = removal(section, grid)
    stats = dict(frac)
    for arm, _ in ARMS:
        stats[arm, "ari"], stats[arm, "agree"] = sm.loc[arm, "ari"], sm.loc[arm, "agree"]
    no = df[df.arm == "01_no_qc"]
    allsq = pd.Index(sorted(set(no.square)))
    draw(section, df, rem, stats, allsq, out)

def run_section(section, grid, out):
    cells = {d: os.path.join(grid, f"{section}_d{d}") for d in (100, 12)}
    raw = {d: ds94().read_bins(c) for d, c in cells.items()}
    names8 = raw[100].obs_names.to_numpy()
    labs, rem, stats = {}, {}, {}
    for arm, _ in ARMS:
        for d in (100, 12):
            f = os.path.join(cells[d], "results", arm, "clean.txt")
            drop = set(open(f).read().split()) if os.path.exists(f) else set()
            keep8 = np.array([b not in drop for b in raw[d].obs_names])
            sub = ds94().pool(raw[d][keep8].copy(), 32)
            labs[arm, d] = ds94().cluster(sub, seed=0)
            rem[arm, d] = removed_fraction(raw[d].obs_names.to_numpy(), keep8)
            stats[arm, d] = 1 - keep8.mean()
            print(f"{section} {arm} d{d}: {100 * stats[arm, d]:.1f}% of 8 µm bins removed, "
                  f"{labs[arm, d].nunique()} clusters on {len(labs[arm, d])} squares", flush=True)
    ref = labs["01_no_qc", 100]
    rows = []
    for arm, _ in ARMS:
        a, b = labs[arm, 100], labs[arm, 12]
        sh = a.index.intersection(b.index)
        ari = adjusted_rand_score(a[sh], b[sh])
        b_on_a = match(b, a)
        agree = (b_on_a[sh] == a[sh])
        for d in (100, 12):
            col = match(labs[arm, d], ref, f"{arm[:2]}d{d}")
            for sq, lab in col.items():
                rows.append((arm, d, sq, lab, float(rem[arm, d].get(sq, 0.0)),
                             bool(agree.get(sq, False)) if sq in sh else None))
        stats[arm, "ari"], stats[arm, "agree"] = ari, agree.mean()
        print(f"  {arm}: ARI {ari:.3f}, {100 * agree.mean():.1f}% of {len(sh)} shared squares agree", flush=True)
    df = pd.DataFrame(rows, columns=["arm", "depth", "square", "label", "removed", "agree"])
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, f"{section}_labels.csv.gz"), index=False)
    allsq = pd.Index(sorted(set(labs["01_no_qc", 100].index) | set(labs["01_no_qc", 12].index)))
    draw(section, df, rem, stats, allsq, out)
    return stats

def palette(labels):
    base = [plt.get_cmap(n)(i) for n in ("tab20", "tab20b", "tab20c") for i in range(20)]
    ref_l = sorted([l for l in labels if not str(l).startswith("new")], key=lambda x: int(x))
    new_l = sorted([l for l in labels if str(l).startswith("new")])
    return {l: base[k % len(base)] for k, l in enumerate(ref_l + new_l)}

SHORT = ROW
W = 7.2
INK, MUTED, GRID = "#1f2a24", "#5a6660", "#e3e8e5"

def draw(section, df, rem, stats, allsq, out):
    r, c = rc_of(allsq)
    R, C = r.max() + 1, c.max() + 1
    pos = {sq: (i, j) for sq, i, j in zip(allsq, r, c)}
    pal = palette(df["label"].unique())
    nrow, asp = len(ARMS), R / C
    fig = plt.figure(figsize=(W, (4 * asp * nrow + 0.6) * W / 18.5))
    outer = fig.add_gridspec(1, 2, width_ratios=[3, 1.25], wspace=0.34,
                             left=0.115, right=0.985, top=0.97, bottom=0.07)
    left = outer[0].subgridspec(nrow, 3, wspace=0.06, hspace=0.20)
    right = outer[1].subgridspec(2, 1, hspace=0.45)
    ax = np.array([[fig.add_subplot(left[k, j]) for j in range(3)] for k in range(nrow)])
    blank = np.full((R, C, 4), 1.0)

    def img(values, colour):
        im = blank.copy()
        for sq in allsq:
            i, j = pos[sq]; im[i, j] = (0.80, 0.80, 0.80, 1)
        for sq, v in values.items():
            if sq in pos:
                i, j = pos[sq]; im[i, j] = colour(v)
        return im

    cmap = plt.get_cmap("Purples")
    for k, (arm, name) in enumerate(ARMS):
        for col, d in enumerate((100, 12)):
            s = df[(df.arm == arm) & (df.depth == d)].set_index("square")["label"]
            ax[k, col].imshow(img(s, lambda v: pal[v]), interpolation="nearest")
            ax[k, col].set_title(f"{d} UMIs per bin, {s.nunique()} clusters", fontsize=5.6, color=INK,
                                 pad=2.5)
        ax[k, 2].imshow(img(rem[arm, 12], lambda v: cmap(0.12 + 0.88 * v)), interpolation="nearest")
        ax[k, 2].set_title("Removed at 12 UMIs per bin", fontsize=5.6, color=INK, pad=2.5)
        ax[k, 0].set_ylabel(ROW[arm], fontsize=6.5, fontweight="bold", rotation=0, ha="right", va="center",
                            labelpad=5, color=INK, linespacing=1.2)
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
    cax = ax[-1, 2].inset_axes([0.05, -0.10, 0.9, 0.04])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 100), cmap=ListedColormap(
        cmap(np.linspace(0.12, 1, 256)))), cax=cax, orientation="horizontal")
    cb.set_ticks([0, 50, 100]); cb.set_ticklabels(["0%", "50%", "100%"])
    cb.set_label("share of the 32 µm square removed", fontsize=5.4, color=MUTED, labelpad=1.5)
    cb.ax.tick_params(labelsize=5.2, colors=MUTED, length=1.5, pad=1); cb.outline.set_visible(False)

    names = [SHORT[a] for a, _ in ARMS]
    y = np.arange(nrow)[::-1]
    r100 = np.array([100 * stats[a, 100] for a, _ in ARMS]); r12 = np.array([100 * stats[a, 12] for a, _ in ARMS])
    ari = np.array([stats[a, "ari"] for a, _ in ARMS])
    c100, c12 = cmap(0.45), cmap(0.92)

    def frame(a, xlabel):
        a.set_xlabel(xlabel, fontsize=5.8, color=MUTED, labelpad=2)
        a.set_yticks(y); a.set_yticklabels(names, fontsize=6, color=INK, linespacing=1.1)
        a.set_ylim(-0.6, nrow - 0.4)
        a.tick_params(axis="x", labelsize=5.4, colors=MUTED, length=0, pad=1.5)
        a.tick_params(axis="y", length=0, pad=2)
        a.grid(axis="x", color=GRID, lw=0.5); a.set_axisbelow(True)
        for sp in ("top", "right", "left"):
            a.spines[sp].set_visible(False)
        a.spines["bottom"].set_color(GRID)

    a1 = fig.add_subplot(right[0])
    for yy, u, v in zip(y, r100, r12):
        a1.plot([u, v], [yy, yy], color="#b9b3d6", lw=1.2, zorder=1, solid_capstyle="round")
        a1.text(max(u, v) + 2.5, yy, f"{u:.1f}% → {v:.1f}%", va="center", fontsize=5.6, color=INK)
    a1.scatter(r100, y, s=14, color=c100, edgecolor="white", linewidth=0.8, zorder=2, label="100 UMIs per bin")
    a1.scatter(r12, y, s=14, color=c12, edgecolor="white", linewidth=0.8, zorder=3, label="12 UMIs per bin")
    a1.set_xlim(0, max(60.0, r12.max() * 1.55))
    frame(a1, "% of 8 µm bins removed")
    a1.legend(frameon=False, fontsize=5.6, loc="upper right", bbox_to_anchor=(1.0, -0.16), ncol=2,
              handletextpad=0.2, columnspacing=0.8, labelcolor=INK)

    a2 = fig.add_subplot(right[1])
    a2.barh(y, ari, height=0.5, color="#5b6b74", zorder=2)
    a2.axvline(ari[0], color=MUTED, lw=0.6, ls=(0, (4, 3)), zorder=3)
    for yy, v in zip(y, ari):
        if v >= 0.08:
            a2.text(v - 0.008, yy, f"{v:.3f}", va="center", ha="right", fontsize=5.6, color="white", fontweight="bold")
        else:
            a2.text(max(v, 0) + 0.008, yy, f"{v:.3f}", va="center", fontsize=5.6, color=INK)
    a2.set_xlim(0, max(0.5, ari.max() * 1.35))
    frame(a2, "ARI between 100 and 12 UMIs per bin")

    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{section}.{ext}"), dpi=600)
    plt.close(fig)
    print(f"  wrote {section}.png/.pdf", flush=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="/work/benchmark_grid")
    ap.add_argument("--sections", required=True)
    ap.add_argument("--out", default="/work/benchmark_grid/figures/demo_downstream")
    ap.add_argument("--replot", action="store_true", help="redraw from saved labels and summary.csv")
    a = ap.parse_args()
    if a.replot:
        for s in a.sections.split(","):
            replot(s, a.grid, a.out)
        return
    summ = []
    for s in a.sections.split(","):
        st = run_section(s, a.grid, a.out)
        for arm, _ in ARMS:
            summ.append(dict(section=s, arm=arm, removed_100=st[arm, 100], removed_12=st[arm, 12],
                             ari=st[arm, "ari"], agree=st[arm, "agree"]))
    pd.DataFrame(summ).to_csv(os.path.join(a.out, "summary.csv"), index=False)
    print("done", flush=True)

if __name__ == "__main__":
    main()

import os

PREFIX = "Visium_HD_"
PROBE, POLYA = "probe", "3prime"
CHEM_LABEL = {PROBE: "Probe-based", POLYA: "3′ polyA"}
CHEMS = (PROBE, POLYA)

SECTIONS = {
    "Human_Colon_Cancer":                   ("human", PROBE, 6.5),
    "Human_Kidney_FFPE":                    ("human", PROBE, 6.5),
    "Human_Lung_Cancer":                    ("human", PROBE, 6.5),
    "Human_Pancreas":                       ("human", PROBE, 6.5),
    "Human_Prostate_Cancer_FFPE":           ("human", PROBE, 6.5),
    "Mouse_Brain":                          ("mouse", PROBE, 6.5),
    "Mouse_Brain_Fresh_Frozen":             ("mouse", PROBE, 6.5),
    "Mouse_Embryo":                         ("mouse", PROBE, 6.5),
    "Mouse_Kidney":                         ("mouse", PROBE, 6.5),
    "Mouse_Small_Intestine":                ("mouse", PROBE, 6.5),
    "Human_Breast_Cancer_Fresh_Frozen":     ("human", PROBE, 6.5),
    "Human_Lung_Cancer_Fixed_Frozen":       ("human", PROBE, 6.5),
    "Human_Tonsil_Fresh_Frozen":            ("human", PROBE, 6.5),
    "Human_Tonsil_Fresh_Frozen_IF":         ("human", PROBE, 6.5),
    "Human_Lymph_Node_FFPE":                ("human", PROBE, 6.5),
    "6p5mm_Human_Heart":                    ("human", PROBE, 6.5),
    "Mouse_Brain_Fixed_Frozen":             ("mouse", PROBE, 6.5),
    "Mouse_Lung_Fresh_Frozen":              ("mouse", PROBE, 6.5),
    "Human_Breast_Cancer_FFPE":             ("human", PROBE, 6.5),
    "Human_Ovarian_Cancer_FF":              ("human", PROBE, 6.5),
    "Human_Colon_Cancer_P1":                ("human", PROBE, 6.5),
    "Human_Colon_Cancer_P5":                ("human", PROBE, 6.5),
    "Human_Colon_Normal_P3":                ("human", PROBE, 6.5),
    "Human_Colon_Normal_P5":                ("human", PROBE, 6.5),
    "SixCancer_Bladder_Cancer":             ("human", PROBE, 6.5),
    "SixCancer_Colorectal_Cancer":          ("human", PROBE, 6.5),
    "SixCancer_Lung_Cancer":                ("human", PROBE, 6.5),
    "SixCancer_Lymphoma":                   ("human", PROBE, 6.5),
    "SixCancer_Ovarian_Cancer":             ("human", PROBE, 6.5),
    "11mm_Human_Breast_Cancer":             ("human", PROBE, 11.0),
    "11mm_Human_Colon_Cancer_Fresh_Frozen": ("human", PROBE, 11.0),
    "11mm_Human_Pancreas":                  ("human", PROBE, 11.0),
    "11mm_Mouse_Embryo":                    ("mouse", PROBE, 11.0),
    "3prime_Human_Lymph_Node":              ("human", POLYA, 6.5),
    "3prime_Human_Ovarian_Cancer":          ("human", POLYA, 6.5),
    "3prime_Human_Tonsil":                  ("human", POLYA, 6.5),
    "3prime_Mouse_Brain":                   ("mouse", POLYA, 6.5),
    "3prime_Mouse_Embryo":                  ("mouse", POLYA, 6.5),
    "3prime_Mouse_Kidney":                  ("mouse", POLYA, 6.5),
    "3prime_Mouse_Ovary":                   ("mouse", POLYA, 6.5),
}
DEPTHS = (100, 50, 25, 12)

CORE_CONDITIONS = ("clean", "dropout01_32um", "dropout01_64um", "dropout01_256um",
                   "mito_32um", "mito_64um", "mito_256um")

SPOTSWEEPER = {8: "24_spotsweeper_8um_k48", 16: "25_spotsweeper_16um_k48",
               32: "26_spotsweeper_32um_k48", 64: "27_spotsweeper_64um_k48"}
SPOTSWEEPER_PACKAGE_DEFAULT_8UM = "03_spotsweeper_8um"
WIDENED = {16: "44_spotsweeper_8um_k192", 32: "45_spotsweeper_8um_k768",
           64: "46_spotsweeper_8um_k3072"}
MAD = {8: "02_mad_8um", 16: "08_mad_16um", 32: "09_mad_32um", 64: "10_mad_64um"}

def short(name):
    return name[len(PREFIX):] if name.startswith(PREFIX) else name

def chemistry(name):
    return SECTIONS[short(name)][1]

def species(name):
    return SECTIONS[short(name)][0]

def cells(grid, depths=DEPTHS, need=("scored", "scores.csv")):
    for sec in sorted(SECTIONS):
        for d in depths:
            path = os.path.join(grid, f"{PREFIX}{sec}_d{d}")
            if os.path.exists(os.path.join(path, *need)):
                yield sec, d, path

MITO_NO_CONTRAST = set()

def scores(cell):
    import pandas as pd
    f = pd.read_csv(os.path.join(cell, "scored", "scores.csv"))
    kind = f["kind"].astype(str)
    stale = kind.str.startswith("mitorel")
    core = os.path.join(cell, "injection_params_mito_relative_core.json")
    if os.path.exists(core) and os.path.getmtime(os.path.join(cell, "scored", "scores.csv")) < os.path.getmtime(core):
        stale |= kind.str.startswith("mito")
    return f[~stale].reset_index(drop=True)

def big512_runs(root):
    p = os.path.join(root, "data", "big512_runs.txt")
    return {ln.strip() for ln in open(p) if ln.strip()} if os.path.exists(p) else set()

def method_failures(grid, arm, big=()):
    import glob
    out = {}
    for d in glob.glob(os.path.join(grid, f"{PREFIX}*_d*")):
        tag = os.path.basename(d)
        try:
            dep = int(tag.rsplit("_d", 1)[1])
        except ValueError:
            continue
        if not os.path.isdir(d):
            continue
        conds = [os.path.basename(c.rstrip("/")) for c in glob.glob(os.path.join(d, "*/"))
                 if os.path.isdir(os.path.join(c, "square_008um"))]
        rdir = os.path.join(d, "results", arm)
        have = set(os.listdir(rdir)) if os.path.isdir(rdir) else set()
        f, n = out.get(dep, (0, 0))
        for c in conds:
            if tag in big and c.endswith("_512um"):
                continue
            n += 1
            log = os.path.join(d, "results", "logs", f"{arm}__{c}.log")
            if f"{c}.txt" not in have and os.path.exists(log) and "ERROR" in open(log, errors="replace").read():
                f += 1
        out[dep] = (f, n)
    return dict(sorted(out.items(), reverse=True))

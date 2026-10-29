import glob
import os

WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARMS = ("01_no_qc,02_mad_8um,08_mad_16um,09_mad_32um,10_mad_64um,03_spotsweeper_8um,"
        "24_spotsweeper_8um_k48,25_spotsweeper_16um_k48,26_spotsweeper_32um_k48,"
        "27_spotsweeper_64um_k48,30_miqc_8um,31_miqc_64um,32_fixedumi_8um,33_fixedumi_64um,"
        "19_spatialartifacts_8um,20_spatialartifacts_16um,34_mad_64um_percovered,"
        "35_mad_64um_excludepartial").split(",")
DEPTHS = {"100", "50", "25", "12"}

try:
    BIG_LIST = {l.strip() for l in open(f"{WK}/data/big512_runs.txt") if l.strip()}
except OSError:
    BIG_LIST = set()
ARMS512 = {"01_no_qc", "02_mad_8um", "08_mad_16um", "09_mad_32um", "10_mad_64um",
           "24_spotsweeper_8um_k48", "25_spotsweeper_16um_k48", "26_spotsweeper_32um_k48",
           "27_spotsweeper_64um_k48"}

RESULT_ERRORS = ("Log-likelihood: NaN",)

def attempted(cell, arm, cond):
    if os.path.exists(f"{cell}/results/{arm}/{cond}.txt"):
        return True
    log = f"{cell}/results/logs/{arm}__{cond}.log"
    if not os.path.exists(log):
        return False
    text = open(log).read()
    return "ERROR" in text and any(e in text for e in RESULT_ERRORS)

for cell in sorted(glob.glob(f"{WK}/benchmark_grid/Visium_HD_*_d*")):
    name, _, depth = os.path.basename(cell).rpartition("_d")
    if depth not in DEPTHS or not os.path.isdir(f"{cell}/results"):
        continue
    conds = [c for c in os.listdir(cell)
             if os.path.isdir(f"{cell}/{c}") and c not in ("results", "scored", "figures")]
    big = os.path.basename(cell) in BIG_LIST
    need = [a for a in ARMS if not all(attempted(cell, a, c) for c in conds
                                       if not (big and c.endswith("_512um") and a not in ARMS512))]
    sc, mk = f"{cell}/scored/scores.csv", f"{cell}/scored/marker.csv"
    stale = (not os.path.exists(sc)
             or any(os.path.getmtime(f) > os.path.getmtime(sc)
                    for f in glob.glob(f"{cell}/results/*/*.txt"))
             or not os.path.exists(mk)
             or os.path.getmtime(sc) > os.path.getmtime(mk))
    if need or stale:
        print(f"{name}\t{depth}\t{','.join(need) if need else '-'}")

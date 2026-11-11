#!/usr/bin/env bash
set -uo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
CACHE=${CACHE:-$WK/cache}
IMG=${IMG:-vhd-qc-benchmark:0.3}
P=${P:-6}
GB=${GB:-60}
LOW_GB=${LOW_GB:-40}
HIGH_GB=${HIGH_GB:-100}
NAME=vhdqc_sa_mad_region
JOBS=$WK/logs/sa_mad_region_jobs.txt
. "$WK/scripts/00_slot.sh"
cd "$WK"
mkdir -p benchmark_saswitch "$CACHE"; chmod 777 "$CACHE"

python3 - "$JOBS" <<'EOF'
import gzip, os, sys
WK = os.environ.get("VHDQC_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{WK}/scripts"); import panel as P
GROUPS = (["clean", "dropout01edge_32um", "dropout01edge_64um", "dropout01edge_256um",
           "mitoedge_32um", "mitoedge_64um", "mitoedge_256um"],
          ["dropout01_32um", "dropout01_64um", "dropout01_256um", "mito_32um", "mito_64um", "mito_256um"])
jobs = []
for sec, dep, cell in P.cells(f"{WK}/benchmark_grid", (100, 12)):
    run = os.path.basename(cell.rstrip("/")); name = run.rsplit("_d", 1)[0]
    if not os.path.isdir(f"{WK}/public/{name}/binned_outputs/square_008um/raw_feature_bc_matrix"):
        print("no raw:", name, file=sys.stderr); continue
    b = 8 if dep == 100 else 16
    out = f"benchmark_saswitch/{run}/47_spatialartifacts_mad_{b}um"
    with gzip.open(f"{cell}/clean/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz", "rt") as fh:
        nb = sum(1 for _ in fh)
    for g, conds in enumerate(GROUPS):
        owed = [c for c in conds if os.path.isdir(f"{cell}/{c}") and not os.path.exists(f"{WK}/{out}/{c}.txt")]
        if not owed:
            continue
        parts = [[c] for c in owed] if (b == 8 and nb > 1_000_000) else [owed]
        for part in parts:
            line = (f"/work/public/{name}/binned_outputs/square_008um /work/benchmark_grid/{run} {b} "
                    f"/work/{out} {','.join(part)}")
            jobs.append((g, 0 if b == 8 else 1, -nb, line))
jobs.sort()
open(sys.argv[1], "w").write("".join(j[-1] + "\n" for j in jobs))
n = [sum(len(j[-1].split()[-1].split(",")) for j in jobs if j[0] == g) for g in (0, 1)]
print(f"{len(jobs)} jobs; condition files owed: edge plan {n[0]}, interior {n[1]}")
EOF
[ -n "${DRY:-}" ] && exit 0

echo "### 154 start $(date -Is), $(wc -l < "$JOBS") jobs, $P workers"
wait_for_slot
docker run --rm --name "$NAME" -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -v "$WK":/work -v "$CACHE":/cache "$IMG" bash -c '
  run() { src=$1; cell=$2; bin=$3; out=$4; conds=$5
    echo "  job start: $(basename $cell) $bin $conds $(date -Is)"
    SOURCE_8UM=$src SA_CACHE=/cache Rscript /work/scripts/153_sa_mad_region.R "$cell" "$bin" "$out" "$conds" 2>&1 |
      grep -E "check1|Error|rror in|Killed" | sed "s/^/  /"
    echo "  job done: $(basename $cell) $bin $conds $(date -Is)"; }
  export -f run
  xargs -P '"$P"' -L 1 bash -c "run \$0 \$1 \$2 \$3 \$4" < /work/logs/sa_mad_region_jobs.txt' &
pid=$!
echo "$(date +%s) $GB $pid" >> "$WK/logs/widen_d100_commit"

guard() {
  local paused="" a
  while kill -0 "$pid" 2>/dev/null; do
    a=$(free -g | awk '/^Mem:/ { print $7 }')
    if [ -z "$paused" ] && [ "$a" -lt "$LOW_GB" ]; then
      docker exec "$NAME" pkill -STOP -x xargs && paused=1 && echo "### dispatch paused, available ${a} GB $(date -Is)"
    elif [ -n "$paused" ] && [ "$a" -gt "$HIGH_GB" ]; then
      docker exec "$NAME" pkill -CONT -x xargs && paused="" && echo "### dispatch resumed, available ${a} GB $(date -Is)"
    fi
    sleep 120
  done
}
guard &
wait "$pid"; rc=$?
echo "### 154 end $(date -Is), exit $rc, condition files $(ls benchmark_saswitch/*/47_spatialartifacts_mad_*/*.txt 2>/dev/null | wc -l)"

set -euo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
. "$WK/scripts/00_slot.sh"
NAME=$1; D=$2; WORKERS=${3:-4}
IMG=${IMG:-vhd-qc-benchmark}
TAG=${NAME}_d${D}; B=/work/benchmark_grid/$TAG
SRC=$WK/public/$NAME/binned_outputs
[ -d "$SRC/square_008um" ] || { echo "no $SRC/square_008um: fetch $NAME first" >&2; exit 1; }
mkdir -p "$WK/benchmark_grid/$TAG"
run() { wait_for_slot; docker run --rm -e OMP_NUM_THREADS=1 -v "$WK":/work -v "$SRC":/src:ro "$IMG" "$@"; }
PY=/opt/pyenv/bin/python
INJ=(--src /src --out "$B" --target-umi-per-bin "$D" --keep-p 0.01 --mito-odds-ratio 10
     --mito-baseline native-mean)
ARMS=01_no_qc,02_mad_8um,08_mad_16um,09_mad_32um,10_mad_64um,03_spotsweeper_8um,\
24_spotsweeper_8um_k48,25_spotsweeper_16um_k48,26_spotsweeper_32um_k48,27_spotsweeper_64um_k48,\
30_miqc_8um,31_miqc_64um,32_fixedumi_8um,33_fixedumi_64um,19_spatialartifacts_8um,\
20_spatialartifacts_16um,34_mad_64um_percovered,35_mad_64um_excludepartial
CORE=clean,dropout01_32um,dropout01_64um,dropout01_256um,mito_32um,mito_64um,mito_256um
EDGE=dropout01edge_32um,dropout01edge_64um,dropout01edge_256um,mitoedge_32um,mitoedge_64um,mitoedge_256um
ARMS512=01_no_qc,02_mad_8um,08_mad_16um,09_mad_32um,10_mad_64um,41_mad_128um,24_spotsweeper_8um_k48,\
25_spotsweeper_16um_k48,26_spotsweeper_32um_k48,27_spotsweeper_64um_k48,28_spotsweeper_128um_k48
WIDE=44_spotsweeper_8um_k192,45_spotsweeper_8um_k768,46_spotsweeper_8um_k3072

run $PY /work/scripts/10_inject_artifacts.py "${INJ[@]}" --diameters 32 64 256 \
    --kinds mito dropout --kind-label dropout01 --with-clean
run $PY /work/scripts/10_inject_artifacts.py "${INJ[@]}" --diameters 32 64 256 \
    --kinds dropout mito --kind-label dropout=dropout01 --location edge --record-suffix _edge
BIG=0; grep -qx "$TAG" "$WK/data/big512_runs.txt" && BIG=1
if [ "$BIG" = 1 ]; then
  run $PY /work/scripts/10_inject_artifacts.py "${INJ[@]}" --diameters 512 --n-regions 8 \
      --kinds dropout mito --kind-label dropout=dropout01 --seed 1 --record-suffix _big512
  run $PY /work/scripts/75_add_coarse_bin.py --benchmark "$B" --bin-um 128
fi

R12() { run Rscript /work/scripts/12_run_arms_k48all.R --benchmark "$B" --workers "$WORKERS" "$@"; }
R12 --arms "$ARMS" --conds "$CORE,$EDGE"
if [ "$BIG" = 1 ]; then
  R12 --arms "$ARMS512" --conds dropout01_512um,mito_512um
  R12 --arms 28_spotsweeper_128um_k48,41_mad_128um --conds clean,dropout01_256um,mito_256um
fi
grep -qx "$TAG" "$WK/data/widened_runs.txt" && R12 --arms "$WIDE" --conds "$CORE"

run $PY /work/scripts/11_score_benchmark.py --benchmark "$B" --results "$B/results" --out "$B/scored"
run $PY /work/scripts/50_marker_validation.py --benchmark "$B" --name "$TAG" --out "$B/scored/marker.csv"
echo "### $TAG done $(date -Is)"

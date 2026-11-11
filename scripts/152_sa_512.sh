#!/usr/bin/env bash
set -uo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
. "$WK/scripts/00_slot.sh"
IMG=${IMG:-vhd-qc-benchmark:0.3}
CONDS="dropout01_512um mito_512um"
LANES=${LANES:-2}
RESERVE_GB=${RESERVE_GB:-80}
SETTLE_S=${SETTLE_S:-60}
Q=$WK/logs/sa_512_queue.txt
LOCK=$WK/logs/widen_d100.lock
COMMIT=$WK/logs/widen_d100_commit
COMMIT_S=${COMMIT_S:-3600}
cd "$WK"

arm_of() { case "${1##*_d}" in 100) echo 19_spatialartifacts_8um ;; 12) echo 20_spatialartifacts_16um ;; esac; }
bins() { zcat "benchmark_grid/$1/clean/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz" | wc -l; }
workers() { if [ "$1" -gt 1000000 ]; then echo 1; else echo 2; fi; }
est_gb() { awk -v b="$1" -v w="$2" 'BEGIN { printf "%d", w * (6 + 40 * b / 1720912) + 1 }'; }
avail_gb() { free -g | awk '/^Mem:/ { print $7 }'; }
committed_gb() {
  local now t g p s=0
  now=$(date +%s)
  if [ -r "$COMMIT" ]; then
    while read -r t g p; do
      [ $((now - t)) -lt "$COMMIT_S" ] && kill -0 "$p" 2>/dev/null && s=$((s + g))
    done < "$COMMIT"
  fi
  echo "$s"
}
fits() { [ $(( $(avail_gb) - $(committed_gb) )) -ge $(($1 + RESERVE_GB)) ]; }
owed() {
  local t=$1 a=$2 c out=""
  for c in $CONDS; do
    [ -d "benchmark_grid/$t/$c/square_008um" ] || continue
    [ -f "benchmark_grid/$t/results/$a/$c.txt" ] && continue
    grep -qs ERROR "benchmark_grid/$t/results/logs/${a}__$c.log" && continue
    out="$out,$c"
  done
  echo "${out#,}"
}

plan() {
  for d in benchmark_grid/*_d100/ benchmark_grid/*_d12/; do
    t=$(basename "$d"); a=$(arm_of "$t")
    [ -n "$(owed "$t" "$a")" ] && echo "$(bins "$t") $t"
  done | sort -k1,1nr
}

if [ -n "${DRY:-}" ]; then
  plan | while read -r b t; do
    a=$(arm_of "$t"); w=$(workers "$b")
    echo "$t $b bins $a w$w ~$(est_gb "$b" "$w")GB [$(owed "$t" "$a")]"
  done
  exit 0
fi

plan | awk '{ print $2 " " $1 }' > "$Q"
echo "### 152 start: $(wc -l < "$Q") runs, lanes $LANES, reserve ${RESERVE_GB} GB $(date -Is)"

next_run() {
  exec 9> "$Q.lock"; flock 9
  local first; first=$(head -1 "$Q"); [ -n "$first" ] && sed -i 1d "$Q"
  flock -u 9; echo "$first"
}

start_container() {
  local log=$1 need=$2 said=""; shift 2
  while :; do
    if fits "$need"; then
      exec 8> "$LOCK"; flock 8
      fits "$need" && break
      flock -u 8
    fi
    [ -z "$said" ] && echo "### waiting for memory: need ${need}+${RESERVE_GB} GB, available $(avail_gb) GB, committed $(committed_gb) GB $(date -Is)"
    said=1; sleep 120
  done
  wait_for_slot
  docker run --rm -e OMP_NUM_THREADS=1 -e OPENBLAS_NUM_THREADS=1 -e MKL_NUM_THREADS=1 "$@" >> "$log" 2>&1 &
  local pid=$!
  echo "$(date +%s) $need $pid" >> "$COMMIT"
  sleep "$SETTLE_S"; flock -u 8
  wait "$pid"
}

lane() {
  local n=$1 job t b a c w need src log
  while job=$(next_run) && [ -n "$job" ]; do
    read -r t b <<< "$job"
    a=$(arm_of "$t"); c=$(owed "$t" "$a"); [ -z "$c" ] && continue
    src=$WK/public/${t%_d*}/binned_outputs; log=logs/sa_512_$t.log
    w=$(workers "$b"); need=$(est_gb "$b" "$w")
    echo "### lane $n: call $t $a w$w ~${need}GB [$c] start $(date -Is)"
    if start_container "$log" "$need" -v "$WK":/work -v "$src":/src:ro "$IMG" \
         Rscript /work/scripts/12_run_arms_k48all.R --benchmark "/work/benchmark_grid/$t" \
         --workers "$w" --arms "$a" --conds "$c"; then
      echo "### lane $n: call $t $a done $(date -Is)"
    else
      echo "!!!! lane $n: $t $a failed $(date -Is), see $log"
    fi
    ( flock 7
      start_container "$log" 8 -v "$WK":/work -v "$src":/src:ro "$IMG" \
        /opt/pyenv/bin/python /work/scripts/11_score_benchmark.py --benchmark "/work/benchmark_grid/$t" \
        --results "/work/benchmark_grid/$t/results" --out "/work/benchmark_grid/$t/scored"
    ) 7> "logs/widen_edge_rescore_$t.lock"
    echo "### lane $n: $t done, rescored $(date -Is)"
  done
  echo "### lane $n: queue empty $(date -Is)"
}

for n in $(seq 1 "$LANES"); do lane "$n" & sleep 5; done
wait
left=0
for d in benchmark_grid/*_d100/ benchmark_grid/*_d12/; do t=$(basename "$d"); [ -n "$(owed "$t" "$(arm_of "$t")")" ] && left=$((left + 1)); done
echo "### 152 end: $left runs still owed $(date -Is)"

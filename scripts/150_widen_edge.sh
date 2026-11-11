#!/usr/bin/env bash
set -uo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
. "$WK/scripts/00_slot.sh"
IMG=${IMG:-vhd-qc-benchmark:0.3}
ARMS=${ARMS:-"44_spotsweeper_8um_k192 45_spotsweeper_8um_k768 46_spotsweeper_8um_k3072"}
CONDS="dropout01edge_32um dropout01edge_64um dropout01edge_256um mitoedge_32um mitoedge_64um mitoedge_256um"
LANES=${LANES:-5}
RESERVE_GB=${RESERVE_GB:-80}
SETTLE_S=${SETTLE_S:-300}
Q=$WK/logs/widen_edge_queue.txt
LOCK=$WK/logs/widen_d100.lock
COMMIT=$WK/logs/widen_d100_commit
MEMLOG=$WK/logs/widen_edge_mem.log
PARAMS=$WK/logs/widen_edge_params
FACTOR=1.9
COMMIT_S=${COMMIT_S:-3600}
cd "$WK"

bins() { zcat "benchmark_grid/$1/clean/square_008um/filtered_feature_bc_matrix/barcodes.tsv.gz" | wc -l; }
k_of() { echo "${1##*_k}"; }
workers() {
  local k; k=$(k_of "$1")
  if [ "$k" -ge 3072 ]; then
    if [ "$2" -gt 1000000 ]; then echo 1; elif [ "$2" -gt 500000 ]; then echo 2; else echo 4; fi
  else echo 4; fi
}
params() { [ -r "$PARAMS" ] && . "$PARAMS"; return 0; }
est_gb() {
  params
  awk -v b="$2" -v k="$(k_of "$1")" -v w="$3" -v f="$FACTOR" 'BEGIN { printf "%d", w * (b * k * 12 * f / 1e9 + 4) + 1 }'
}
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

prio() {
  case "$1:${2##*_d}" in
    44_*:100|45_*:100) echo 1 ;;
    44_*:12|45_*:12) echo 2 ;;
    *) echo 3 ;;
  esac
}

plan() {
  for d in benchmark_grid/*_d100/ benchmark_grid/*_d12/; do
    t=$(basename "$d"); b=""
    for a in $ARMS; do
      [ -z "$(owed "$t" "$a")" ] && continue
      [ -z "$b" ] && b=$(bins "$t")
      echo "$(prio "$a" "$t") $b $t $a"
    done
  done | sort -k1,1n -k2,2n
}

if [ -n "${DRY:-}" ]; then
  plan | while read -r p b t a; do
    w=$(workers "$a" "$b")
    echo "p$p $t $b bins $(k_of "$a")nb w$w ~$(est_gb "$a" "$b" "$w")GB [$(owed "$t" "$a" | tr ',' '\n' | wc -l)]"
  done
  exit 0
fi

if [ -z "${JOIN:-}" ]; then
  plan | awk '{ print $3 " " $2 " " $4 }' > "$Q"
  echo "### 150 start: $(wc -l < "$Q") run-arm pairs, arms $ARMS, image $IMG, lanes $LANES, reserve ${RESERVE_GB} GB $(date -Is)"
  monitor() {
    while :; do
      docker ps --format '{{.ID}} {{.Image}}' < /dev/null 2>/dev/null | awk '$2 ~ /^(vhd-qc-benchmark|vhdqc-wknd)/ { print $1 }' |
        while read -r id; do
          echo "$(date -Is) $id $(docker stats --no-stream --format '{{.MemUsage}}' "$id" < /dev/null 2>/dev/null) avail=$(avail_gb)G"
        done >> "$MEMLOG"
      sleep 120
    done
  }
  monitor & MON=$!
  trap 'kill $MON 2>/dev/null' EXIT
else
  MON=
  echo "### 150 join: $(wc -l < "$Q") run-arm pairs left in the queue, $LANES more lanes $(date -Is)"
fi

next_run() {
  exec 9> "$Q.lock"; flock 9
  local first; first=$(head -1 "$Q"); [ -n "$first" ] && sed -i 1d "$Q"
  flock -u 9; echo "$first"
}

start_container() {
  local log=$1 need=$2 settle=$3 said=""; shift 3
  while :; do
    params
    if fits "$need"; then
      exec 8> "$LOCK"; flock 8; params
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
  sleep "${settle:-$SETTLE_S}"; flock -u 8
  wait "$pid"
}

lane() {
  local n=$1 job t b a c w need src log
  while job=$(next_run) && [ -n "$job" ]; do
    read -r t b a <<< "$job"
    src=$WK/public/${t%_d*}/binned_outputs; log=logs/widen_edge_$t.log
    c=$(owed "$t" "$a"); [ -z "$c" ] && continue
    w=$(workers "$a" "$b"); need=$(est_gb "$a" "$b" "$w")
    echo "### lane $n: call $t $a w$w ~${need}GB start $(date -Is)"
    if start_container "$log" "$need" "" -v "$WK":/work -v "$src":/src:ro "$IMG" \
         Rscript /work/scripts/12_run_arms_k48all.R --benchmark "/work/benchmark_grid/$t" \
         --workers "$w" --arms "$a" --conds "$c"; then
      echo "### lane $n: call $t $a done $(date -Is)"
    else
      echo "!!!! lane $n: $t $a failed $(date -Is), see $log"
    fi
    echo "### lane $n: call $t rescore start $(date -Is)"
    ( flock 7
      start_container "$log" 8 30 -v "$WK":/work -v "$src":/src:ro "$IMG" \
        /opt/pyenv/bin/python /work/scripts/11_score_benchmark.py --benchmark "/work/benchmark_grid/$t" \
        --results "/work/benchmark_grid/$t/results" --out "/work/benchmark_grid/$t/scored"
    ) 7> "logs/widen_edge_rescore_$t.lock"
    echo "### lane $n: $t $a done, rescored $(date -Is)"
  done
  echo "### lane $n: queue empty $(date -Is)"
}

for n in $(seq 1 "$LANES"); do lane "${JOIN:+j}$n" & sleep 5; done
wait $(jobs -p | grep -v "^${MON:-x}$")
left=0
for d in benchmark_grid/*_d100/ benchmark_grid/*_d12/; do for a in $ARMS; do [ -n "$(owed "$(basename "$d")" "$a")" ] && left=$((left + 1)); done; done
echo "### 150 end: $left run-arm pairs still owed $(date -Is)"

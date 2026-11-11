#!/usr/bin/env bash
set -uo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
cd "$WK"
. "$WK/scripts/00_slot.sh"
export PRIORITY=1
IMG=vhd-qc-benchmark:0.4
TMP=benchmark_grid/figures/downstream_guideline
LANES=${LANES:-3}
NEED_GB=${NEED_GB:-30}
RESERVE_GB=${RESERVE_GB:-80}
COMMIT=$WK/logs/widen_d100_commit
COMMIT_S=3600
ONLY=${ONLY:-}
mkdir -p "$TMP"

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

SECS=()
if [ -n "$ONLY" ]; then SECS=("$ONLY")
else
  for d in benchmark_grid/*_d100; do
    s=$(basename "$d" _d100)
    [ -d "benchmark_grid/${s}_d12" ] && SECS+=("$s")
  done
fi
echo "### 158 start: ${#SECS[@]} sections, $LANES lanes, image $IMG, need ${NEED_GB}+${RESERVE_GB} GB $(date -Is)"

lane() {
  local n=$1 i s pid said
  for ((i = n; i < ${#SECS[@]}; i += LANES)); do
    s=${SECS[$i]}
    [ -s "$TMP/$s.csv" ] && { echo "### lane $n: $s skip (done)"; continue; }
    said=""
    until fits "$NEED_GB"; do
      [ -z "$said" ] && echo "### lane $n: waiting for memory, available $(avail_gb) GB, committed $(committed_gb) GB $(date -Is)"
      said=1; sleep 120
    done
    wait_for_slot
    docker run --rm --name "vhdqc_dsguide_$n" -e OMP_NUM_THREADS=6 -v "$WK":/work "$IMG" \
      /opt/pyenv/bin/python /work/scripts/158_downstream_guideline.py --section "$s" \
      --out "/work/$TMP/$s.csv" > "logs/downstream_guideline_$s.log" 2>&1 &
    pid=$!
    echo "$(date +%s) $NEED_GB $pid" >> "$COMMIT"
    echo "### lane $n: $s start $(date -Is)"
    if wait "$pid"; then echo "### lane $n: $s ok $(date -Is)"
    else echo "!!!! lane $n: $s FAILED $(date -Is), see logs/downstream_guideline_$s.log"; fi
  done
}
for ((n = 0; n < LANES; n++)); do lane "$n" & sleep 20; done
wait
echo "### per-section files: $(ls "$TMP"/*.csv 2>/dev/null | wc -l) of ${#SECS[@]}"
echo "### 158 end $(date -Is)"

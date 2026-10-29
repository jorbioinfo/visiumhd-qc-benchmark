WK=${WK:-${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}
CAP=${CAP:-5}
CAP_FILE=${CAP_FILE:-$WK/logs/slot_cap}
slot_cap() {
  local n e
  if [ -r "$CAP_FILE" ] && read -r n e < "$CAP_FILE" && [[ "${n:-}" =~ ^[0-9]+$ ]] && [ -n "${e:-}" ] &&
     [ "$(date +%s)" -lt "$(date -d "$e" +%s 2>/dev/null || echo 0)" ]; then echo "$n"; else echo "$CAP"; fi
}
SLOT_DIR=${SLOT_DIR:-$WK/logs/slot_priority}
CLAIM_S=${CLAIM_S:-120}
my_containers() {
  local ids n=0
  ids=$(docker ps -q 2>/dev/null) || return 0
  [ -z "$ids" ] && { echo 0; return 0; }
  while read -r id; do
    case "$(docker inspect -f '{{.Config.Image}}' "$id" < /dev/null 2>/dev/null)" in
      vhd-qc-benchmark*|vhdqc-wknd*) n=$((n+1)) ;;
    esac
  done <<< "$ids"
  echo "$n"
}
priority_waiting() {
  local f m n=0 now
  now=$(date +%s)
  for f in "$SLOT_DIR"/*; do
    [ -e "$f" ] || continue
    case "${f##*/}" in claim.*) continue ;; esac
    m=$(stat -L -c %Y "$f" 2>/dev/null) || continue
    [ $((now - m)) -le 90 ] && n=$((n+1))
  done
  echo "$n"
}
claims_fresh() {
  local f m n=0 now
  now=$(date +%s)
  for f in "$SLOT_DIR"/claim.*; do
    [ -e "$f" ] || continue
    m=$(stat -c %Y "$f" 2>/dev/null) || continue
    if [ $((now - m)) -le "$CLAIM_S" ]; then n=$((n+1)); else rm -f "$f"; fi
  done
  echo "$n"
}
wait_for_slot() {
  local n p=0 c tok="" got
  mkdir -p "$SLOT_DIR"
  [ "${PRIORITY:-0}" = 1 ] && tok=$SLOT_DIR/prio.$BASHPID
  while :; do
    [ -n "$tok" ] && touch "$tok"
    got=$(
      exec 9> "$SLOT_DIR/.lock"
      flock -w 300 9 || { echo "lock timeout"; exit; }
      n=$(my_containers); c=$(claims_fresh); cap=$(slot_cap)
      [ -z "$tok" ] && p=$(priority_waiting)
      if [ $((n + c + p)) -lt "$cap" ]; then
        touch "$SLOT_DIR/claim.$BASHPID.$(date +%s%N)"; echo ok
      else
        echo "$n/$cap in use, $c just claimed, $p held for priority"
      fi
    )
    if [ "$got" = ok ]; then [ -n "$tok" ] && rm -f "$tok"; return 0; fi
    echo "   slot: ${got}${tok:+ (priority)}, waiting $(date -Is)"
    sleep 60
  done
}

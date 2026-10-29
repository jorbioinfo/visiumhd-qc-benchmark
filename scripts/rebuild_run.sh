set -uo pipefail
WK=${VHDQC_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}
. "$WK/scripts/00_slot.sh"
RUN=$1; shift
NAME=${RUN%_d*}
IMG=${IMG:-vhd-qc-benchmark}
wait_for_slot
docker run --rm -e OMP_NUM_THREADS=1 -e ONLY_GIVEN="${ONLY_GIVEN:-}" -v "$WK":/work -v "$WK/public/$NAME/binned_outputs":/src:ro \
  "$IMG" /opt/pyenv/bin/python /work/scripts/136_rebuild_check.py "$RUN" "$@"

#!/usr/bin/env bash
# Runs query-index with several VEOs on the same index and queries, and summarizes the times.
# Usage: [TIMEOUT=600] [QUERY_INDEX=...] scripts/bench.sh <index> <queries> [limit=1000] [out_dir=bench-out] [veo ...]
# Output: <out_dir>/<veo>.out (query;results;nanoseconds per query) and <out_dir>/summary.csv
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
bin="${QUERY_INDEX:-$here/../build/query-index}"
index="$1"; queries="$2"; limit="${3:-1000}"; out="${4:-bench-out}"
shift $(( $# < 4 ? $# : 4 ))
veos=("$@")
if [ ${#veos[@]} -eq 0 ]; then
  veos=(global adaptive global-int3 global-guided3
        adaptive-int1 adaptive-int2 adaptive-int3 adaptive-int4
        adaptive-guided1 adaptive-guided2 adaptive-guided3 adaptive-guided4)
fi
mkdir -p "$out"
for v in "${veos[@]}"; do
  echo "running $v" >&2
  "$bin" "$index" "$queries" "$limit" "$v" "${TIMEOUT:-600}" | grep ';' > "$out/$v.out"
done
python3 "$here/summarize.py" "$out" "${veos[@]}" | tee "$out/summary.csv"

#!/usr/bin/env bash
# Builds the compressed (rrr) and plain indexes of a dataset and benchmarks all VEOs on both.
# Usage: scripts/run-all.sh <file.dat> <queries> [limit=1000] [out_dir=results]
# Queries with a variable predicate are not supported by this index and are filtered out first.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; bin="$here/../build"
dat="$1"; queries="$2"; limit="${3:-1000}"; out="${4:-results}"
mkdir -p "$out"
awk '{ok=1; n=split($0, t, " "); for (i=2; i<=n; i+=4) if (substr(t[i],1,1)=="?") ok=0; if (ok) print}' "$queries" > "$out/queries.fixedp.txt"
echo "$(wc -l < "$out/queries.fixedp.txt") of $(wc -l < "$queries") queries have constant predicates" | tee "$out/info.txt"
[ -f "$dat.ring" ] || "$bin/build-index" "$dat" | tee -a "$out/info.txt"
[ -f "$dat.ring-plain" ] || "$bin/build-index-plain" "$dat" | tee -a "$out/info.txt"
ls -l "$dat.ring" "$dat.ring-plain" | tee -a "$out/info.txt"
veos=(global adaptive global-int3 global-guided3 adaptive-int1 adaptive-int2 adaptive-int3
      adaptive-guided1 adaptive-guided2 adaptive-guided3)
QUERY_INDEX="$bin/query-index" "$here/bench.sh" "$dat.ring" "$out/queries.fixedp.txt" "$limit" "$out/rrr" "${veos[@]}"
QUERY_INDEX="$bin/query-index-plain" "$here/bench.sh" "$dat.ring-plain" "$out/queries.fixedp.txt" "$limit" "$out/plain" "${veos[@]}"
# estimator cost (time spent computing weights) for the adaptive variants
for idx in rrr plain; do
  suffix=$([ $idx = plain ] && echo "-plain" || echo "")
  ix=$([ $idx = plain ] && echo "$dat.ring-plain" || echo "$dat.ring")
  for v in adaptive adaptive-int1 adaptive-int3 adaptive-guided3; do
    "$bin/query-index$suffix-stats" "$ix" "$out/queries.fixedp.txt" "$limit" "$v" | grep ';' > "$out/$idx/stats.$v.out"
    awk -F';' -v v="$v" -v ix="$idx" '{t+=$3; c+=$4; e+=$5} END{printf "%s,%s,%.1f,%d,%.1f\n", ix, v, t/1e6, c, e/1e6}' "$out/$idx/stats.$v.out"
  done
done | (echo "index,veo,total_ms,weight_computations,weight_ms"; cat) | tee "$out/estimator-cost.csv"

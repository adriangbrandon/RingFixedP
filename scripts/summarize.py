#!/usr/bin/env python3
# Summarizes the outputs of bench.sh: average/median time per VEO, speedup over the first VEO, and a check that
# all VEOs report the same number of results per query.
import sys, statistics
out, veos = sys.argv[1], sys.argv[2:]
data = {}
for v in veos:
    rows = [l.strip().split(';') for l in open('%s/%s.out' % (out, v)) if l.strip()]
    data[v] = {int(q): (int(n), int(t) / 1e6) for q, n, t in rows}
base = veos[0]
print('veo,queries,avg_ms,median_ms,p90_ms,max_ms,speedup_avg_vs_%s,median_ratio_vs_%s,wins_vs_%s,count_mismatches' % (base, base, base))
for v in veos:
    ts = [t for _, t in data[v].values()]
    ts_sorted = sorted(ts)
    mism = sum(1 for q in data[v] if q in data[base] and data[v][q][0] != data[base][q][0])
    ratios = [data[base][q][1] / max(data[v][q][1], 1e-6) for q in data[v] if q in data[base]]
    wins = sum(1 for q in data[v] if q in data[base] and data[v][q][1] < data[base][q][1])
    print('%s,%d,%.3f,%.3f,%.3f,%.3f,%.2f,%.2f,%d,%d' % (
        v, len(ts), statistics.mean(ts), statistics.median(ts), ts_sorted[int(0.9 * (len(ts) - 1))], max(ts),
        statistics.mean(data[base][q][1] for q in data[v]) / statistics.mean(ts), statistics.median(ratios), wins, mism))

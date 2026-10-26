import glob,os,sys,statistics,math
V=['adaptive','adaptive-int1','adaptive-int2','adaptive-int3','adaptive-guided1','adaptive-guided2','adaptive-guided3']
def load(res,ix):
    runs={}
    for v in V:
        r={}
        for l in open(f'{res}/2d-diego-{ix}/{v}.out'):
            p=l.strip().split(';')
            if len(p)>=3 and p[0].isdigit(): r[int(p[0])]=int(p[2])/1e6
        runs[v]=r
    return runs
for label,res in [(a,b) for a,b in zip(sys.argv[1::2],sys.argv[2::2])]:
  for ix in ['rrr','plain']:
    R=load(res,ix); qs=sorted(R[V[0]])
    best={q:min(R[v][q] for v in V[1:]) for q in qs}
    print(f'\n{label} {ix}: variant, avg_ms, total_s, p99_ms, max_s, geomean_vs_best, bad(>2x best & >+100ms), sum_excess_over_best_s')
    for v in V:
        r=R[v]; ts=sorted(r.values())
        g=math.exp(statistics.mean(math.log(max(r[q],1e-3)/max(best[q],1e-3)) for q in qs))
        bad=sum(1 for q in qs if r[q]>2*best[q] and r[q]-best[q]>100)
        exc=sum(r[q]-best[q] for q in qs)/1e3
        print(f'{v:18s} {statistics.mean(ts):8.1f} {sum(ts)/1e3:7.0f} {ts[int(.99*(len(ts)-1))]:9.0f} {ts[-1]/1e3:7.1f} {g:6.3f} {bad:4d} {exc:7.0f}')

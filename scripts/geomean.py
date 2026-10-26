import glob,os,sys,statistics,math
res=sys.argv[1]; runs={}
for f in glob.glob(res+'/2d-*/*.out')+glob.glob(res+'/3d/*.out'):
    if os.path.basename(f).startswith('stats.'): continue
    n=os.path.relpath(f,res)[:-4].replace('/',':'); r={}
    for l in open(f):
        p=l.strip().split(';')
        if len(p)>=3 and p[0].isdigit(): r[int(p[0])]=int(p[2])/1e6
    if r: runs[n]=r
base=sys.argv[2]; b=runs[base]
print('vs',base)
for n in sys.argv[3:]:
    if n not in runs: print(n,'missing'); continue
    r=runs[n]; qs=[q for q in r if q in b]
    g=math.exp(statistics.mean(math.log(max(b[q],1e-3)/max(r[q],1e-3)) for q in qs))
    print('%-45s total %8.1f s  geomean %.2f'%(n,sum(r[q] for q in qs)/1e3,g))

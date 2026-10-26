# Synthetic fixed-predicate graph + queries + brute-force answers.
import random, sys, itertools, collections
random.seed(int(sys.argv[2]) if len(sys.argv) > 2 else 7)
out = sys.argv[1]
U, P, N = 1000, 6, 8000
# skewed node popularity
w = [1.0 / (i + 1) ** 0.8 for i in range(U)]
nodes = list(range(1, U + 1))
T = set()
while len(T) < N:
    p = random.choices(range(1, P + 1), weights=[8, 5, 3, 2, 1, 1])[0]
    s = random.choices(nodes, weights=w)[0] if random.random() < .5 else random.randint(1, U)
    o = random.choices(nodes, weights=w)[0] if random.random() < .5 else random.randint(1, U)
    T.add((s, p, o))
T = sorted(T)
with open(out + '.dat', 'w') as f:
    f.write('\n'.join('%d %d %d' % t for t in T))   # no trailing newline (build-index reads until eof)
byp = collections.defaultdict(list)
for s, p, o in T: byp[p].append((s, o))
shapes = [
    "?x {a} ?y . ?y {b} ?z",
    "?x {a} ?y . ?y {b} ?z . ?z {c} ?x",
    "?x {a} ?y . ?x {b} ?z . ?x {c} ?w",
    "?x {a} ?y . ?y {b} ?z . ?z {c} ?w",
    "?x {a} ?y . ?y {b} ?x",
    "?x {a} ?y . ?z {b} ?y . ?x {c} ?z",
    "?x {a} ?y . ?y {b} ?z . ?x {c} ?w . ?w {a} ?z",
    "{n} {a} ?y . ?y {b} ?z",
    "?x {a} ?y . ?y {b} {n}",
    "?x {a} {n} . ?x {b} ?y . ?y {c} ?z",
]
Q = []
for sh in shapes:
    for _ in range(4):
        a, b, c = (random.randint(1, P) for _ in range(3))
        n = random.choice(byp[a])[0] if sh.startswith('{n}') else random.choice(T)[2]
        Q.append(sh.format(a=a, b=b, c=c, n=n))
with open(out + '.queries', 'w') as f:
    f.write('\n'.join(Q) + '\n')
# brute force: DFS over patterns, choosing at each step the pattern with most bound terms
bys = collections.defaultdict(lambda: collections.defaultdict(list))
byo = collections.defaultdict(lambda: collections.defaultdict(list))
pairs = {p: set(v) for p, v in byp.items()}
for p, v in byp.items():
    for s, o in v:
        bys[p][s].append(o); byo[p][o].append(s)
def solve(q):
    pats = [t.split() for t in q.split(' . ')]
    def val(t, mu):
        if not t.startswith('?'): return int(t)
        return mu.get(t)
    def rec(rem, mu):
        if not rem: return 1
        # pick pattern with most bound terms
        best = max(rem, key=lambda i: (val(pats[i][0], mu) is not None) + (val(pats[i][2], mu) is not None))
        s, p, o = pats[best]; p = int(p)
        rest = [i for i in rem if i != best]
        vs, vo = val(s, mu), val(o, mu)
        tot = 0
        if vs is not None and vo is not None:
            return rec(rest, mu) if (vs, vo) in pairs[p] else 0
        if vs is not None:
            cands = [(vs, x) for x in bys[p][vs]]
        elif vo is not None:
            cands = [(x, vo) for x in byo[p][vo]]
        else:
            cands = byp[p]
        for cs, co in cands:
            m = dict(mu)
            if s.startswith('?'):
                if m.setdefault(s, cs) != cs: continue
            if o.startswith('?'):
                if m.setdefault(o, co) != co: continue
            tot += rec(rest, m)
        return tot
    return rec(list(range(len(pats))), {})
with open(out + '.expected', 'w') as f:
    for i, q in enumerate(Q):
        f.write('%d;%d\n' % (i, solve(q)))

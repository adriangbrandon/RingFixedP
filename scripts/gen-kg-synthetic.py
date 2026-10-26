# Larger synthetic knowledge-graph-like dataset with fixed-predicate queries (shapes as in WGPB).
# Nodes belong to classes; each predicate links a source class to a target class; popularity is Zipfian.
import sys, numpy as np
out, N, U, P, NC = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), 30
rng = np.random.default_rng(11)
clustered = len(sys.argv) > 5 and sys.argv[5] == 'clustered'
if clustered:   # identifiers of each class are contiguous (ids correlated with types, as often in real KGs)
    node_class = np.sort(rng.integers(0, NC, U))
else:
    node_class = rng.integers(0, NC, U)
members = [np.flatnonzero(node_class == c) + 1 for c in range(NC)]
for m in members: rng.shuffle(m)   # popularity is unrelated to the identifier inside a class
pf = 1.0 / np.arange(1, P + 1) ** 1.1; pf /= pf.sum()
src, tgt = rng.integers(0, NC, P), rng.integers(0, NC, P)
def zipf_pick(m, k, a):
    # power-law-like popularity: index = floor(len * u^a), larger a = more skew
    r = (len(m) * rng.random(k) ** a).astype(np.int64)
    return m[np.minimum(r, len(m) - 1)]
counts = rng.multinomial(int(N * 1.15), pf)
S, Pp, O = [], [], []
for p in range(P):
    k = counts[p]
    if k == 0: continue
    a_s, a_o = rng.uniform(1.5, 4.0), rng.uniform(2.0, 6.0)
    S.append(zipf_pick(members[src[p]], k, a_s)); O.append(zipf_pick(members[tgt[p]], k, a_o)); Pp.append(np.full(k, p + 1))
T = np.stack([np.concatenate(S), np.concatenate(Pp), np.concatenate(O)], 1).astype(np.int64)
T = np.unique(T, axis=0)
if len(T) > N: T = T[np.sort(rng.choice(len(T), N, replace=False))]
np.savetxt(out + '.dat', T, fmt='%d', delimiter=' ')
with open(out + '.dat', 'rb+') as f:          # build-index reads until eof: drop the trailing newline
    f.seek(-1, 2); f.truncate()
print('triples', len(T))
# adjacency
bs = T[np.lexsort((T[:, 2], T[:, 1], T[:, 0]))]; s_keys = bs[:, 0]
bo = T[np.lexsort((T[:, 0], T[:, 1], T[:, 2]))]; o_keys = bo[:, 2]
pairs = set(map(tuple, T[:, [0, 2]].tolist())) if len(T) < 3_000_000 else None
def out_edges(x):
    a, b = np.searchsorted(s_keys, x), np.searchsorted(s_keys, x, 'right'); return bs[a:b]
def in_edges(x):
    a, b = np.searchsorted(o_keys, x), np.searchsorted(o_keys, x, 'right'); return bo[a:b]
def step(x, fwd):
    e = out_edges(x) if fwd else in_edges(x)
    if len(e) == 0: return None
    t = e[rng.integers(len(e))]
    return (t[1], t[2]) if fwd else (t[1], t[0])
def edge_between(x, y):
    e = out_edges(x); e = e[e[:, 2] == y]
    return None if len(e) == 0 else e[rng.integers(len(e))][1]
def fmt(t, consts):
    return ' '.join(str(consts[v]) if v in consts else '?' + v for v in t)
Q = []
def emit(pats, emb):
    consts = {}
    if rng.random() < 0.25:                    # one constant node
        v = list(emb)[rng.integers(len(emb))]; consts[v] = emb[v]
    Q.append(' . '.join('%s %d %s' % (('%d' % consts[s]) if s in consts else '?' + s, p, ('%d' % consts[o]) if o in consts else '?' + o) for s, p, o in pats))
def rand_triple():
    return T[rng.integers(len(T))]
tries = 0
target = {'path2': 50, 'path3': 50, 'path4': 40, 'star3': 40, 'star4': 30, 'tri': 40, 'sq': 40}
got = dict.fromkeys(target, 0)
while any(got[k] < target[k] for k in target) and tries < 400000:
    tries += 1
    kind = [k for k in target if got[k] < target[k]][rng.integers(sum(got[k] < target[k] for k in target))]
    s, p, o = rand_triple()
    if kind.startswith('path'):
        L = int(kind[-1]); pats = [('x0', p, 'x1')]; emb = {'x0': s, 'x1': o}; cur = o; ok = True
        for i in range(1, L):
            fwd = rng.random() < 0.7
            r = step(cur, fwd)
            if r is None: ok = False; break
            q, nxt = r; a, b = 'x%d' % i, 'x%d' % (i + 1)
            pats.append((a, q, b) if fwd else (b, q, a)); emb[b] = nxt; cur = nxt
        if ok and len(set(emb.values())) == len(emb): emit(pats, emb); got[kind] += 1
    elif kind.startswith('star'):
        L = int(kind[-1]); e = out_edges(s)
        if len(e) < L: continue
        sel = e[rng.choice(len(e), L, replace=False)]
        if len(set(sel[:, 1])) < 2: continue
        pats = [('x', t[1], 'y%d' % i) for i, t in enumerate(sel)]; emb = {'x': s}; emb.update({'y%d' % i: t[2] for i, t in enumerate(sel)})
        emit(pats, emb); got[kind] += 1
    elif kind == 'tri':
        r = step(o, True)
        if r is None: continue
        q, z = r
        c = edge_between(s, z)
        if c is None or len({s, o, z}) < 3: continue
        emit([('x', p, 'y'), ('y', q, 'z'), ('x', c, 'z')], {'x': s, 'y': o, 'z': z}); got[kind] += 1
    else:
        r = step(o, True)
        if r is None: continue
        q, z = r
        e = in_edges(z); e = e[e[:, 0] != o]
        if len(e) == 0: continue
        t = e[rng.integers(len(e))]; w = t[0]
        c = edge_between(s, w)
        if c is None or len({s, o, z, w}) < 4: continue
        emit([('x', p, 'y'), ('y', q, 'z'), ('x', c, 'w'), ('w', t[1], 'z')], {'x': s, 'y': o, 'z': z, 'w': w}); got[kind] += 1
print('queries', got, 'tries', tries)
rng.shuffle(Q)
with open(out + '.queries', 'w') as f: f.write('\n'.join(Q) + '\n')

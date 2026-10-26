# Adaptive VEOs and refined estimators for the Ring with fixed predicates

This report explains what was added to RingFixedP (the Ring with fixed predicates of Arroyuelo et al., TODS 2024,
Section 4.4) and how it compares on Wikidata with the original order and with the 3D Ring / IRing of
[ring-basic-p](https://github.com/adriangbrandon/ring-basic-p). The summaries of the runs are in
`doc/results/wikidata/`.

## 1. Summary

- RingFixedP now has the two improvements of Section 6 of "New Compressed Indices for Multijoins on Graph Databases"
  (Information Systems 2025, the IRing paper): an **adaptive variable elimination order (VEO)**, where only the next
  variable is chosen, after each binding; and the **refined estimator of the size of the intersections** (Eq. 5),
  which sums the minimum number of triples of the patterns over a fixed number of parts of the identifier space,
  instead of taking the minimum size of the intervals (Eq. 1).
- Triple patterns may now have **constant subjects and objects** (before, only predicates could be constant).
- Lonely variables (those that appear in a single pattern) are enumerated lazily, so the result limit applies to
  them; before, queries like `?x P31 Q5` materialized all the values and took up to 80 s.
- **Default: `adaptive-guided3`** (adaptive VEO, Eq. 5 over the 8 quantiles of the smallest interval). On Wikidata
  it is 1.6x (filtered) and 3.9x (full) faster in total time than the original order on the compressed index, and
  2.2x / 5.8x on the plain-bitvector index. It is the best or within run-to-run noise of the best of all the
  fixed-part variants in the four settings.
- The sdsl-lite fork of adriangbrandon lacks the `rrr_vector<15>` specialization, which makes every
  `rrr_vector<15>` (this index and the 3D c-ring) about 2x slower; the one-line fix is in the README.

## 2. What was implemented

### 2.1 Background
RingFixedP stores, for each predicate p, the (subject, object) pairs of p as a wavelet matrix over n_p points:
positions follow the object order and values follow the subject order. A bound object restricts the positions, a
bound subject restricts the values. Leapfrog Triejoin (LTJ) runs on top, with one iterator per triple pattern.

### 2.2 VEOs (`include/veo.hpp`)
- `global`: the original order of RingFixedP (`gao_size_v2`): variables sorted by weight once, before LTJ, related
  variables first, lonely variables last.
- `adaptive`: after each binding, the weights of the unbound variables related to the bound one are recomputed and
  the next variable is the one of minimum weight (IRing Sec. 6.1). Weights are recomputed lazily and restored on
  backtracking, and nothing is estimated when a single candidate remains.

### 2.3 Estimators
- Eq. 1 (`global`, `adaptive`): minimum size of the intervals of the patterns that contain the variable.
- Eq. 5 (`*-int<k>`, `*-guided<k>`, k = 1..4): the identifier space [0, U] is cut into 2^k parts and the weight is
  the sum over the parts of the minimum, over the patterns, of the number of triples of the pattern whose value of
  the variable falls in the part. In the 2D ring every predicate has its own grid, so the parts are defined on
  identifiers (not on the nodes of a shared wavelet matrix, as in the IRing) and the counts are 2D range counts
  (`sub_ring::count_S_lt`, `count_O_lt`, O(log n_p) each).
  - `int<k>` uses 2^k equal parts: the dyadic ranges that the IRing obtains from the first k levels of its wavelet
    matrices.
  - `guided<k>` uses the 2^k quantiles of the smallest interval, found with `kth_S` (one descent) and `kth_O`
    (binary search over counts). Since the smallest interval bounds the estimate, cutting it into parts of equal
    size gives every part the same chance to lower it.

### 2.4 Other changes
- Constant subjects and objects (`ltj_iterator`).
- Lazy enumeration of lonely variables with `leap` (see above).
- Index construction in O(n_p log n_p) per predicate instead of O(P * U); the indexes are byte-identical.
- `build-index-plain` / `query-index-plain`: the same index with plain bitvectors (`RING_PLAIN_BV`).
- `cmake -DSDSL_PREFIX=...` selects the sdsl-lite installation; it compiles with the darroyue and adriangbrandon
  forks.
- `query-index-stats`: also reports the number of weight computations and the time spent in them.
- `tests/test-estimators.cpp`: brute-force checks of the counting and quantile primitives.

## 3. How it was evaluated

- Machine: 2x Intel Xeon Silver 4110 (32 cores), 735 GB RAM, Debian, g++ 14.2, darroyue/sdsl-lite (and
  adriangbrandon/sdsl-lite with the fix for the 3D indexes). One query at a time.
- Data: Wikidata as in the IRing paper. *Filtered* (~80M triples) with the 850 queries of
  `Queries-wikidata-benchmark.txt`; *full* (~1000M triples) with the 1148 queries of `Queries-bgps-limit1000.txt`
  that have constant predicates (out of 1295).
- At most 1000 results per query, timeout 600 s, wall-clock time per query.
- 3D: ring-basic-p with its plain `.ring` and compressed `.c-ring` (rrr, rebuilt with the fixed sdsl), with the
  fixed VEO of the Ring (`query-index-fixed`), the adaptive VEO with Eq. 1 (`query-index`) and the IRing
  (`query-index-veo-intersection`, adaptive VEO with Eq. 5 on k = 3 levels).
- Correctness: all the configurations return the same number of results per query as the original order (except
  3D queries that hit the timeout); brute-force counts on synthetic graphs; `test-estimators`.
- Run-to-run noise on unchanged configurations is about 5-15% of the totals.

## 4. Results

Index sizes:

| index | filtered | full |
|---|---|---|
| 2D rrr (`.ring`) | 399 MB | 4.02 GB |
| 2D plain (`.ring-plain`) | 528 MB | 5.94 GB |
| 3D c-ring | 545 MB | 7.00 GB |
| 3D ring | 909 MB | 11.65 GB |

Total time in seconds over all the queries (from `doc/results/wikidata/*-compare.csv`).

| configuration | filtered rrr | filtered plain | full rrr | full plain |
|---|---|---|---|---|
| global (original) | 50.6 | 21.4 | 1407 (1 timeout) | 713 |
| adaptive (Eq. 1) | 45.1 | 20.4 | 1109 (1 timeout) | 606 |
| global-int3 | 41.0 | 17.2 | 920 | 311 |
| global-guided3 | 28.4 | 10.6 | 913 | 286 |
| adaptive-int1 | 44.3 | 17.0 | 355 | 119 |
| adaptive-int2 | 39.4 | 15.7 | 384 | 123 |
| adaptive-int3 | 43.5 | 15.4 | 398 | 113 |
| adaptive-guided1 | 36.2 | 16.8 | 402 | 124 |
| adaptive-guided2 | 32.7 | 12.2 | 396 | 124 |
| **adaptive-guided3** | **31.8** | **9.9** | **362** | **123** |

3D, same queries:

| configuration | filtered | full |
|---|---|---|
| ring, Ring (fixed VEO) | 44.2 | 1278 (1 timeout) |
| ring, adaptive (Eq. 1) | 46.6 | 1031 (1 timeout) |
| ring, IRing | 45.4 | 401 |
| c-ring, Ring (fixed VEO) | 134.7 | 3336 (2 timeouts) |
| c-ring, adaptive (Eq. 1) | 122.3 | 2409 (1 timeout) |
| c-ring, IRing | 116.8 | 1522 |

Per query (geometric mean of the ratios), `adaptive-guided3` on the 2D rrr index is 1.5x (filtered) / 1.25x (full)
faster than the IRing on the 3D ring, at a third of its size, and 4.1x / 3.7x faster than the IRing on the 3D
c-ring, at 57% of its size. On the 2D plain index it is 3.9x / 2.7x faster than the IRing on the 3D ring, at half
its size. The 3D ring, of course, also supports variable predicates.

Main findings:

1. **The combination matters.** On full Wikidata the adaptive VEO alone gains 1.3x and the refined estimator with
   the global VEO 1.5x (rrr), but together they gain 3.9x (rrr) and 5.8x (plain): they avoid the catastrophic
   orders (one query goes from the 600 s timeout to 2 s).
2. **guided beats int.** Quantiles of the smallest interval give better estimates with fewer parts than equal
   parts of [0, U]: on filtered, guided3 is 1.4x (rrr) / 1.6x (plain) faster than int3. On full all the adaptive
   refined variants are within noise of each other (the total is dominated by a few expensive queries whose order
   they all get right). int1 is not better on full beyond noise, and it is clearly worse on filtered.
3. **The estimator has a cost.** With fixed parts, computing the weights takes 12-34% of the time on filtered
   (guided cuts need `kth_O`, a binary search of counts), and about 1% on full, where the queries are more
   expensive. This cost grows with 2^k and limits how large k can usefully be.
4. **sdsl.** Use darroyue/sdsl-lite, or adriangbrandon/sdsl-lite with the `rrr_vector_15.hpp` include; otherwise
   the compressed indexes are ~2x slower.

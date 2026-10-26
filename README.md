# RingFixedP
A versión of the Ring Index for queries with P constant

This version adds the improvements of Section 6 of "New Compressed Indices for Multijoins on Graph Databases"
(Information Systems, 2025): an **adaptive variable elimination order (VEO)**, where the next variable is chosen
after each binding, and the **refined estimator of the size of the intersections** (Eq. 5), computed over a fixed
number of parts. It also supports constant subjects and objects in the triple patterns, and an optional variant
with plain bitvectors. See `doc/REPORT.md` for the details and the experiments on Wikidata.

## Building

    mkdir build && cd build
    cmake .. [-DSDSL_PREFIX=<prefix>]     # sdsl-lite installed in <prefix>/include, <prefix>/lib (default: $HOME)
    make

It compiles with https://github.com/darroyue/sdsl-lite and https://github.com/adriangbrandon/sdsl-lite. With the
latter, add `#include "rrr_vector_15.hpp"` before the final `#endif` of `include/sdsl/rrr_vector.hpp` (as in the
original sdsl-lite) and reinstall it; otherwise `rrr_vector<15>` does not use its specialization and the compressed
index is about 2x slower (and its file format changes, so indexes must be rebuilt after the fix).

## Usage

    ./build-index data.dat                     # -> data.dat.ring        (compressed, rrr bitvectors)
    ./build-index-plain data.dat               # -> data.dat.ring-plain  (plain bitvectors: ~1.4x larger, ~3x faster)
    ./query-index data.dat.ring queries.txt [limit=1000 (0: no limit)] [veo=adaptive-guided3] [timeout=600]
    ./query-index-plain data.dat.ring-plain queries.txt ...

`data.dat` has one triple `s p o` of integer identifiers per line. Each line of `queries.txt` is a basic graph
pattern, with triple patterns separated by ` . ` and variables starting with `?`; the predicates must be constant.
The output has one line `query;results;nanoseconds` per query.

VEOs (`veo` argument):

| veo | order | estimator |
|---|---|---|
| `global` | fixed before running LTJ (original RingFixedP order) | minimum interval size (Eq. 1) |
| `adaptive` | next variable chosen after each binding | Eq. 1 |
| `{global,adaptive}-int<k>`, k = 1..4 | | Eq. 5 over 2^k equal parts of the identifier space |
| `{global,adaptive}-guided<k>`, k = 1..4 | | Eq. 5 over the 2^k quantiles of the smallest interval |

The default, `adaptive-guided3`, is 1.6x (filtered) and 3.9x (full Wikidata) faster in total time than the original
order on the compressed index, and 2.2x / 5.8x on the plain one.

Other executables: `query-index-stats` / `query-index-plain-stats` also print the number of weight computations and
the time spent in them; `test-estimators` checks the counting primitives against brute force.
`scripts/run-all.sh <data.dat> <queries>` builds both indexes and benchmarks all the VEOs
(`scripts/bench.sh`, `scripts/summarize.py`, `scripts/compare-refined.py`, `scripts/geomean.py`).

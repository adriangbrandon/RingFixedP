// Checks count_S_lt, count_O_lt, kth_S and kth_O of the ring with fixed predicates against a brute force over the
// triples. Usage: test-estimators <file.dat> (builds the index in memory)
#include <iostream>
#include <random>
#include <ring_fixed_p.hpp>

int main(int argc, char *argv[]) {
    std::vector<spo_triple> D;
    std::ifstream ifs(argv[1]);
    uint64_t s, p, o;
    while (ifs >> s >> p >> o) D.emplace_back(s, p, o);
    std::vector<spo_triple> T = D;
    ring::ring_fixed_p R(D);
    uint64_t U = 0, P = 0;
    for (auto &t : T) {
        U = std::max({U, (uint64_t) std::get<0>(t), (uint64_t) std::get<2>(t)});
        P = std::max(P, (uint64_t) std::get<1>(t));
    }
    std::mt19937_64 rng(3);
    uint64_t errors = 0, checks = 0;
    for (uint64_t pv = 1; pv <= P; ++pv) {
        std::vector<std::pair<uint64_t, uint64_t>> so;
        for (auto &t : T) if (std::get<1>(t) == pv) so.emplace_back(std::get<0>(t), std::get<2>(t));
        if (so.empty()) continue;
        ring::bwt_interval full = R.down_P(pv);
        for (int it = 0; it < 300; ++it) {
            // subject variable, object bound (interval of C_s) or unbound
            bool bound = it % 3 != 0;
            uint64_t ob = so[rng() % so.size()].second, sb = so[rng() % so.size()].first;
            ring::bwt_interval I_s = bound ? R.down_P_O(full, pv, ob) : full;
            ring::bwt_interval I_o = bound ? R.down_P_S(full, sb) : full;
            uint64_t v = rng() % (U + 2);
            std::vector<uint64_t> subs, objs;
            for (auto &x : so) {
                if (!bound || x.second == ob) subs.push_back(x.first);
                if (!bound || x.first == sb) objs.push_back(x.second);
            }
            std::sort(subs.begin(), subs.end());
            std::sort(objs.begin(), objs.end());
            uint64_t cs = std::lower_bound(subs.begin(), subs.end(), v) - subs.begin();
            uint64_t co = std::lower_bound(objs.begin(), objs.end(), v) - objs.begin();
            checks += 2;
            if (R.count_S_lt(I_s, v) != cs) ++errors;
            if (R.count_O_lt(I_o, v) != co) ++errors;
            uint64_t q = rng() % subs.size();
            checks++;
            if (R.kth_S(I_s, q) != subs[q]) ++errors;
            q = rng() % objs.size();
            checks++;
            if (R.kth_O(I_o, q) != objs[q]) ++errors;
        }
    }
    std::cout << checks << " checks, " << errors << " errors" << std::endl;
    return errors != 0;
}

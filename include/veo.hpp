/*
 * veo.hpp
 *
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// Variable elimination orders (VEOs) for the ring with fixed predicates, following Section 6 of
// "New Compressed Indices for Multijoins on Graph Databases" (Arroyuelo et al.):
//  - veo_global:   the order is fixed before running LTJ (the original gao_size_v2 behaviour).
//  - veo_adaptive: only the next variable is chosen, after each binding (Section 6.1).
// Both take an estimator of the weight of a variable:
//  - estimator_size:      minimum size of the intervals of the triple patterns (Eq. 1).
//  - estimator_intersect: sum over a partition of [0, U] of the minimum number of triples per part (Eq. 5).
//                         In the fixed-predicate ring the parts must be defined on the identifiers, because each
//                         predicate has its own grid. The weight of each part is a count of points in the grid.
//                         The parts are either 2^k equal intervals of [0, U] (what the wavelet matrices of the
//                         ring compute on k levels) or, if guided, the 2^k quantiles of the smallest interval.

#ifndef RING_VEO_HPP
#define RING_VEO_HPP

#include <ring_fixed_p.hpp>
#include <ltj_iterator.hpp>
#include <triple_pattern.hpp>
#include <utils.hpp>
#include <algorithm>
#include <chrono>
#include <list>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace ring {

    namespace veo {

        // Optional instrumentation (compile with -DVEO_STATS): number of weight computations and time spent in them
        struct stats {
            static uint64_t &calls() { static uint64_t c = 0; return c; }
            static uint64_t &ns() { static uint64_t t = 0; return t; }
            static void reset() { calls() = 0; ns() = 0; }
        };

#ifdef VEO_STATS
        struct stats_timer {
            std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
            ~stats_timer() {
                ++stats::calls();
                stats::ns() += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - t0).count();
            }
        };
#define VEO_STATS_TIMER stats_timer _veo_stats_timer;
#else
#define VEO_STATS_TIMER
#endif

        struct estimator_size {

            template<class ring_t, class iter_t, class var_t>
            static uint64_t weight(ring_t *ring, const var_t x, const std::vector<iter_t *> &itrs) {
                VEO_STATS_TIMER
                uint64_t w = -1ULL;
                for (const iter_t *iter : itrs) {
                    w = std::min(w, util::get_size_interval(*iter));
                }
                return w;
            }
        };

        template<uint64_t k = 3, bool guided = false>
        struct estimator_intersect {

            template<class ring_t>
            static inline uint64_t count_lt(ring_t *ring, bool subject, bwt_interval &I, uint64_t v) {
                return subject ? ring->count_S_lt(I, v) : ring->count_O_lt(I, v);
            }

            template<class ring_t, class iter_t, class var_t>
            static uint64_t weight(ring_t *ring, const var_t x, const std::vector<iter_t *> &itrs) {
                if (k == 0) return estimator_size::weight(ring, x, itrs);
                VEO_STATS_TIMER

                //1. Interval of each triple pattern where the values of x are (C_s for subjects, C_o for objects)
                const uint64_t n = itrs.size();
                std::vector<bool> subject(n);
                std::vector<bwt_interval> intervals(n);
                uint64_t guide = 0;
                for (uint64_t i = 0; i < n; ++i) {
                    subject[i] = itrs[i]->is_variable_subject(x);
                    intervals[i] = subject[i] ? itrs[i]->i_s : itrs[i]->i_o;
                    if (intervals[i].size() == 0) return 0;
                    if (intervals[i].size() < intervals[guide].size()) guide = i;
                }

                //2. Cuts c_0 = 0 < c_1 < ... < c_m = U+1; the j-th part is [c_j, c_{j+1}-1]
                const uint64_t U = ring->universe(intervals[0].get_p_value());
                std::vector<uint64_t> cuts;
                if (!guided) {
                    uint64_t bits = 64 - __builtin_clzll(U);
                    uint64_t step = 1ULL << ((bits > k) ? bits - k : 0);
                    for (uint64_t c = 0; c <= U; c += step) cuts.push_back(c);
                } else {
                    cuts.push_back(0);
                    const uint64_t size = intervals[guide].size(), m = 1ULL << k;
                    for (uint64_t j = 1; j < m; ++j) {
                        uint64_t q = j * size / m;
                        if (q == 0) continue;
                        cuts.push_back(subject[guide] ? ring->kth_S(intervals[guide], q)
                                                      : ring->kth_O(intervals[guide], q));
                    }
                    std::sort(cuts.begin(), cuts.end());
                    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
                }
                cuts.push_back(U + 1);

                //3. Sum over the parts of the minimum number of triples in the part
                std::vector<uint64_t> prev(n, 0);
                uint64_t w = 0;
                for (uint64_t j = 1; j < cuts.size(); ++j) {
                    uint64_t min_part = -1ULL;
                    for (uint64_t i = 0; i < n; ++i) {
                        uint64_t cur = count_lt(ring, subject[i], intervals[i], cuts[j]);
                        min_part = std::min(min_part, cur - prev[i]);
                        prev[i] = cur;
                    }
                    w += min_part;
                }
                return w;
            }
        };

        template<class ltj_iterator_t, class estimator_t = estimator_size>
        class veo_global {

        public:
            typedef ltj_iterator_t ltj_iter_type;
            typedef typename ltj_iter_type::var_type var_type;
            typedef typename ltj_iter_type::ring_type ring_type;
            typedef uint64_t size_type;
            typedef std::unordered_map<var_type, std::vector<ltj_iter_type *>> var_to_iterators_type;
            typedef struct {
                var_type name;
                size_type weight;
                size_type n_triples;
                std::unordered_set<var_type> related;
            } info_var_type;
            typedef std::pair<size_type, var_type> pair_type;
            typedef std::priority_queue<pair_type, std::vector<pair_type>, greater<pair_type>> min_heap_type;

        private:
            std::vector<var_type> m_order;
            size_type m_index = 0;

            struct compare_var_info {
                inline bool operator()(const info_var_type &linfo, const info_var_type &rinfo) {
                    if (linfo.n_triples > 1 && rinfo.n_triples == 1) {
                        return true;
                    }
                    if (linfo.n_triples == 1 && rinfo.n_triples > 1) {
                        return false;
                    }
                    return linfo.weight < rinfo.weight;
                }
            };

            static void fill_heap(const var_type var,
                                  std::unordered_map<var_type, size_type> &hash_table,
                                  std::vector<info_var_type> &vec,
                                  std::vector<bool> &checked,
                                  min_heap_type &heap) {
                auto pos_var = hash_table[var];
                for (const auto &e : vec[pos_var].related) {
                    auto pos_rel = hash_table[e];
                    if (!checked[pos_rel] && vec[pos_rel].n_triples > 1) {
                        heap.push({vec[pos_rel].weight, e});
                        checked[pos_rel] = true;
                    }
                }
            }

        public:
            veo_global() = default;

            // Same order as gao::gao_size_v2, with the weights given by estimator_t
            veo_global(const std::vector<triple_pattern> *triple_patterns,
                       const std::vector<ltj_iter_type> *iterators,
                       const var_to_iterators_type *var_iterators,
                       ring_type *r) {

                //1. Filling var_info with data about each variable (in order of appearance)
                std::vector<info_var_type> var_info;
                std::unordered_map<var_type, size_type> hash_table_position;
                auto add_var = [&](const var_type var) {
                    if (hash_table_position.find(var) != hash_table_position.end()) return;
                    const auto &itrs = var_iterators->at(var);
                    info_var_type info;
                    info.name = var;
                    info.n_triples = itrs.size();
                    info.weight = (itrs.size() > 1) ? estimator_t::weight(r, var, itrs)
                                                    : estimator_size::weight(r, var, itrs);
                    var_info.emplace_back(info);
                    hash_table_position.insert({var, var_info.size() - 1});
                };
                for (const triple_pattern &triple : *triple_patterns) {
                    if (triple.s_is_variable()) add_var((var_type) triple.term_s.value);
                    if (triple.o_is_variable()) add_var((var_type) triple.term_o.value);
                    if (triple.s_is_variable() && triple.o_is_variable()) {
                        var_type var_s = (var_type) triple.term_s.value, var_o = (var_type) triple.term_o.value;
                        var_info[hash_table_position[var_s]].related.insert(var_o);
                        var_info[hash_table_position[var_o]].related.insert(var_s);
                    }
                }

                //2. Sorting variables according to their weights
                std::sort(var_info.begin(), var_info.end(), compare_var_info());
                size_type lonely_start = var_info.size();
                for (size_type i = 0; i < var_info.size(); ++i) {
                    hash_table_position[var_info[i].name] = i;
                    if (var_info[i].n_triples == 1 && i < lonely_start) {
                        lonely_start = i;
                    }
                }

                //3. Choosing the variables: related ones first, lonely ones at the end
                size_type i = 0;
                std::vector<bool> checked(var_info.size(), false);
                m_order.reserve(var_info.size());
                while (i < lonely_start) {
                    if (!checked[i]) {
                        m_order.push_back(var_info[i].name);
                        checked[i] = true;
                        min_heap_type heap;
                        auto var_name = var_info[i].name;
                        fill_heap(var_name, hash_table_position, var_info, checked, heap);
                        while (!heap.empty()) {
                            var_name = heap.top().second;
                            heap.pop();
                            m_order.push_back(var_name);
                            fill_heap(var_name, hash_table_position, var_info, checked, heap);
                        }
                    }
                    ++i;
                }
                while (i < var_info.size()) {
                    m_order.push_back(var_info[i].name);
                    ++i;
                }
            }

            inline var_type next() { return m_order[m_index++]; }

            inline void down() {}

            inline void up() {}

            inline void done() { --m_index; }

            inline size_type size() const { return m_order.size(); }
        };

        // Port of veo::veo_adaptive of ring-basic-p (the ring of the paper)
        template<class ltj_iterator_t, class estimator_t = estimator_size>
        class veo_adaptive {

        public:
            typedef ltj_iterator_t ltj_iter_type;
            typedef typename ltj_iter_type::var_type var_type;
            typedef typename ltj_iter_type::ring_type ring_type;
            typedef uint64_t size_type;
            typedef std::unordered_map<var_type, std::vector<ltj_iter_type *>> var_to_iterators_type;
            typedef struct {
                var_type name;
                size_type weight;
                std::unordered_set<var_type> related;
                bool is_bound;
                bool dirty; //weight not up to date with the current bindings
            } info_var_type;
            typedef struct {
                size_type pos;
                size_type w;
                bool dirty;
            } update_type;
            typedef std::vector<update_type> version_type;

        private:
            const var_to_iterators_type *m_ptr_var_iterators = nullptr;
            ring_type *m_ptr_ring = nullptr;
            std::unordered_map<var_type, size_type> m_hash_table_position;
            std::vector<info_var_type> m_var_info;
            std::vector<var_type> m_lonely;
            std::list<size_type> m_not_bound;
            std::stack<size_type> m_bound;
            std::stack<version_type> m_versions;
            size_type m_index = 0;

            void add_var(const var_type var) {
                if (m_hash_table_position.find(var) != m_hash_table_position.end()) return;
                if (std::find(m_lonely.begin(), m_lonely.end(), var) != m_lonely.end()) return;
                const auto &itrs = m_ptr_var_iterators->at(var);
                if (itrs.size() == 1) {
                    m_lonely.emplace_back(var);
                    return;
                }
                info_var_type info;
                info.name = var;
                info.weight = estimator_t::weight(m_ptr_ring, var, itrs);
                info.is_bound = false;
                info.dirty = false;
                m_var_info.emplace_back(info);
                m_hash_table_position.insert({var, m_var_info.size() - 1});
                m_not_bound.push_back(m_var_info.size() - 1);
            }

        public:
            veo_adaptive() = default;

            veo_adaptive(const std::vector<triple_pattern> *triple_patterns,
                         const std::vector<ltj_iter_type> *iterators,
                         const var_to_iterators_type *var_iterators,
                         ring_type *r) {
                m_ptr_var_iterators = var_iterators;
                m_ptr_ring = r;
                for (const triple_pattern &triple : *triple_patterns) {
                    if (triple.s_is_variable()) add_var((var_type) triple.term_s.value);
                    if (triple.o_is_variable()) add_var((var_type) triple.term_o.value);
                    if (triple.s_is_variable() && triple.o_is_variable()) {
                        auto it_s = m_hash_table_position.find((var_type) triple.term_s.value);
                        auto it_o = m_hash_table_position.find((var_type) triple.term_o.value);
                        if (it_s != m_hash_table_position.end() && it_o != m_hash_table_position.end()) {
                            m_var_info[it_s->second].related.insert(m_var_info[it_o->second].name);
                            m_var_info[it_o->second].related.insert(m_var_info[it_s->second].name);
                        }
                    }
                }
            }

            // The unbound (non-lonely) variable with the smallest weight; then the lonely ones.
            // Weights are recomputed lazily: only the stale ones, and only if there is more than one candidate.
            inline var_type next() {
                if (m_index < m_var_info.size()) {
                    size_type min = -1ULL;
                    auto min_iter = m_not_bound.begin();
                    const bool single = (m_not_bound.size() == 1);
                    for (auto iter = m_not_bound.begin(); iter != m_not_bound.end() && !single; ++iter) {
                        auto &info = m_var_info[*iter];
                        if (info.dirty) {
                            info.weight = estimator_t::weight(m_ptr_ring, info.name, m_ptr_var_iterators->at(info.name));
                            info.dirty = false;
                        }
                        if (min > info.weight) {
                            min = m_var_info[*iter].weight;
                            min_iter = iter;
                        }
                    }
                    size_type min_pos = *min_iter;
                    m_var_info[min_pos].is_bound = true;
                    m_bound.emplace(min_pos);
                    m_not_bound.erase(min_iter);
                    ++m_index;
                    return m_var_info[min_pos].name;
                } else {
                    ++m_index;
                    return m_lonely[m_index - m_var_info.size() - 1];
                }
            }

            // Called after binding the last variable returned by next(): the weights of its related variables become stale
            inline void down() {
                if (m_index - 1 < m_var_info.size()) {
                    version_type version;
                    for (const auto &rel : m_var_info[m_bound.top()].related) {
                        const auto pos = m_hash_table_position[rel];
                        if (!m_var_info[pos].is_bound) {
                            version.emplace_back(update_type{pos, m_var_info[pos].weight, m_var_info[pos].dirty});
                            m_var_info[pos].dirty = true;
                        }
                    }
                    m_versions.emplace(version);
                }
            }

            // Called before unbinding: restores the weights
            inline void up() {
                if (m_index - 1 < m_var_info.size()) {
                    for (const update_type &update : m_versions.top()) {
                        m_var_info[update.pos].weight = update.w;
                        m_var_info[update.pos].dirty = update.dirty;
                    }
                    m_versions.pop();
                }
            }

            // Called when all the values of the last variable returned by next() have been processed
            inline void done() {
                --m_index;
                if (m_index < m_var_info.size()) {
                    auto pos = m_bound.top();
                    m_not_bound.push_back(pos);
                    m_var_info[pos].is_bound = false;
                    m_bound.pop();
                }
            }

            inline size_type size() const { return m_var_info.size() + m_lonely.size(); }
        };
    }
}

#endif //RING_VEO_HPP

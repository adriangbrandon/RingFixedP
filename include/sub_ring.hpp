/*
 * sub_ring.hpp
 * Copyright (C) 2020 Author name removed for double blind evaluation
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

#ifndef SUB_RING
#define SUB_RING

#include <cstdint>
#include <configuration.hpp>
#include <fstream>

class sub_ring {

public:
    typedef uint64_t size_type;
    typedef uint64_t value_type;

private:
    table_column_type m_F_s;
    uint64_t m_n_p;  // number of tuples

    sdsl::sd_vector<> m_B_s;
    sdsl::rank_support_sd<> m_rank_B_s;
    sdsl::select_support_sd<> m_select_B_s;
    /**/
    sdsl::sd_vector<> m_A_s;
    sdsl::rank_support_sd<> m_rank_A_s;
    sdsl::select_support_sd<> m_select_A_s;
    /**/
    sdsl::sd_vector<> m_B_o;
    sdsl::rank_support_sd<> m_rank_B_o;
    sdsl::select_support_sd<> m_select_B_o;
    /**/
    sdsl::sd_vector<> m_A_o;
    sdsl::rank_support_sd<> m_rank_A_o;
    sdsl::select_support_sd<> m_select_A_o;


    void copy(const sub_ring &o) {
        m_F_s = o.m_F_s;
        m_n_p = o.m_n_p;

        m_B_s = o.m_B_s;
        m_rank_B_s = o.m_rank_B_s;
        m_rank_B_s.set_vector(&m_B_s);
        m_select_B_s = o.m_select_B_s;
        m_select_B_s.set_vector(&m_B_s);

        m_A_s = o.m_A_s;
        m_rank_A_s = o.m_rank_A_s;
        m_rank_A_s.set_vector(&m_A_s);
        m_select_A_s = o.m_select_A_s;
        m_select_A_s.set_vector(&m_A_s);

        m_B_o = o.m_B_o;
        m_rank_B_o = o.m_rank_B_o;
        m_rank_B_o.set_vector(&m_B_o);
        m_select_B_o = o.m_select_B_o;
        m_select_B_o.set_vector(&m_B_o);

        m_A_o = o.m_A_o;
        m_rank_A_o = o.m_rank_A_o;
        m_rank_A_o.set_vector(&m_A_o);
        m_select_A_o = o.m_select_A_o;
        m_select_A_o.set_vector(&m_A_o);
    }


public:
    sub_ring() = default;


    // Builds a subring from a vector of pairs
    // Builds a subring from a vector of (subject, object) pairs. Takes O(n_p log n_p) time: the arrays of size U of
    // the first version (O(P*U) time overall) are replaced by the distinct subjects/objects (same index produced).
    sub_ring(std::vector<std::pair<uint64_t, uint64_t>> &D, uint64_t U) {
        cout << D.size() << "  ";
        fflush(stdout);
        m_n_p = D.size();
        // First sort using the first component of the pairs
        sort(D.begin(), D.end(), [](const std::pair<uint64_t, uint64_t> &a,
                                    const std::pair<uint64_t, uint64_t> &b) { return a.first < b.first; });
        // Distinct subjects and the number of triples with a smaller subject (A_s)
        std::vector<uint64_t> subjects, starts_s;
        for (uint64_t i = 0; i < D.size(); i++) {
            if (i == 0 || D[i].first != D[i - 1].first) {
                subjects.push_back(D[i].first);
                starts_s.push_back(i);
            }
        }

        // Then, stable sort using the second component of the pairs: the order of C_s
        stable_sort(D.begin(), D.end(), [](const std::pair<uint64_t, uint64_t> &a,
                                           const std::pair<uint64_t, uint64_t> &b) { return a.second < b.second; });
        // Distinct objects and the number of triples with a smaller object (A_o)
        std::vector<uint64_t> objects, starts_o;
        for (uint64_t i = 0; i < D.size(); i++) {
            if (i == 0 || D[i].second != D[i - 1].second) {
                objects.push_back(D[i].second);
                starts_o.push_back(i);
            }
        }

        // F_s[i] = A_s[C_s[i]] + rank of C_s[i] in C_s[1..i] (LF)
        sdsl::int_vector<> _F_s(D.size() + 2, 0);
        {
            std::vector<uint64_t> seen(subjects.size(), 0);
            for (uint64_t i = 0; i < D.size(); i++) {
                uint64_t j = std::lower_bound(subjects.begin(), subjects.end(), D[i].first) - subjects.begin();
                _F_s[i + 1] = starts_s[j] + (++seen[j]);
            }
        }
        util::bit_compress(_F_s);
        construct_im(m_F_s, _F_s);

        // Bitvectors: B marks the identifiers that appear (plus U+1), A marks the starting positions (plus n)
        auto build_sd = [](const std::vector<uint64_t> &ones, uint64_t last, uint64_t size) {
            sdsl::sd_vector_builder builder(size, ones.size() + 1);
            for (const auto &x : ones) builder.set(x);
            builder.set(last);
            return sdsl::sd_vector<>(builder);
        };
        m_B_s = build_sd(subjects, U + 1, U + 2);
        m_A_s = build_sd(starts_s, D.size(), D.size() + 1);
        m_B_o = build_sd(objects, U + 1, U + 2);
        m_A_o = build_sd(starts_o, D.size(), D.size() + 1);
        util::init_support(m_rank_B_s, &m_B_s);
        util::init_support(m_select_B_s, &m_B_s);
        util::init_support(m_rank_A_s, &m_A_s);
        util::init_support(m_select_A_s, &m_A_s);
        util::init_support(m_rank_B_o, &m_B_o);
        util::init_support(m_select_B_o, &m_B_o);
        util::init_support(m_rank_A_o, &m_A_o);
        util::init_support(m_select_A_o, &m_A_o);

        cout << "-- Subring constructed successfully" << endl;
        fflush(stdout);
    };

    //! Copy constructor
    sub_ring(const sub_ring &o) {
        copy(o);
    }

    //! Move constructor
    sub_ring(sub_ring &&o) {
        *this = std::move(o);
    }

    //! Copy Operator=
    sub_ring &operator=(const sub_ring &o) {
        if (this != &o) {
            copy(o);
        }
        return *this;
    }

    //! Move Operator=
    sub_ring &operator=(sub_ring &&o) {
        if (this != &o) {
            m_F_s = o.m_F_s;
            m_n_p = o.m_n_p;

            m_B_s = std::move(o.m_B_s);
            m_rank_B_s = std::move(o.m_rank_B_s);
            m_rank_B_s.set_vector(&m_B_s);
            m_select_B_s = std::move(o.m_select_B_s);
            m_select_B_s.set_vector(&m_B_s);

            m_A_s = std::move(o.m_A_s);
            m_rank_A_s = std::move(o.m_rank_A_s);
            m_rank_A_s.set_vector(&m_A_s);
            m_select_A_s = std::move(o.m_select_A_s);
            m_select_A_s.set_vector(&m_A_s);

            m_B_o = std::move(o.m_B_o);
            m_rank_B_o = std::move(o.m_rank_B_o);
            m_rank_B_o.set_vector(&m_B_o);
            m_select_B_o = std::move(o.m_select_B_o);
            m_select_B_o.set_vector(&m_B_o);

            m_A_o = std::move(o.m_A_o);
            m_rank_A_o = std::move(o.m_rank_A_o);
            m_rank_A_o.set_vector(&m_A_o);
            m_select_A_o = std::move(o.m_select_A_o);
            m_select_A_o.set_vector(&m_A_o);
        }
        return *this;
    }

    void swap(sub_ring &o) {

        std::swap(m_F_s, o.m_F_s);
        std::swap(m_n_p, o.m_n_p);

        std::swap(m_B_s, o.m_B_s);
        sdsl::util::swap_support(m_rank_B_s, o.m_rank_B_s, &m_B_s, &o.m_B_s);
        sdsl::util::swap_support(m_select_B_s, o.m_select_B_s, &m_B_s, &o.m_B_s);
        std::swap(m_A_s, o.m_A_s);
        sdsl::util::swap_support(m_rank_A_s, o.m_rank_A_s, &m_A_s, &o.m_A_s);
        sdsl::util::swap_support(m_select_A_s, o.m_select_A_s, &m_A_s, &o.m_A_s);

        std::swap(m_B_o, o.m_B_o);
        sdsl::util::swap_support(m_rank_B_o, o.m_rank_B_o, &m_B_o, &o.m_B_o);
        sdsl::util::swap_support(m_select_B_o, o.m_select_B_o, &m_B_o, &o.m_B_o);
        std::swap(m_A_o, o.m_A_o);
        sdsl::util::swap_support(m_rank_A_o, o.m_rank_A_o, &m_A_o, &o.m_A_o);
        sdsl::util::swap_support(m_select_A_o, o.m_select_A_o, &m_A_o, &o.m_A_o);
    }

    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += sdsl::write_member(m_n_p, out, child, "n_p");
        written_bytes += m_F_s.serialize(out, child, "F_s");
        written_bytes += m_B_s.serialize(out, child, "B_s");
        written_bytes += m_rank_B_s.serialize(out, child, "rank_B_s");
        written_bytes += m_select_B_s.serialize(out, child, "select_B_s");
        written_bytes += m_A_s.serialize(out, child, "A_s");
        written_bytes += m_rank_A_s.serialize(out, child, "rank_A_s");
        written_bytes += m_select_A_s.serialize(out, child, "select_A_s");
        written_bytes += m_B_o.serialize(out, child, "B_o");
        written_bytes += m_rank_B_o.serialize(out, child, "rank_B_o");
        written_bytes += m_select_B_o.serialize(out, child, "select_B_o");
        written_bytes += m_A_o.serialize(out, child, "A_o");
        written_bytes += m_rank_A_o.serialize(out, child, "rank_A_o");
        written_bytes += m_select_A_o.serialize(out, child, "select_A_o");
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream &in) {

        sdsl::read_member(m_n_p, in);
        m_F_s.load(in);
        m_B_s.load(in);
        m_rank_B_s.load(in);
        m_rank_B_s.set_vector(&m_B_s);
        m_select_B_s.load(in);
        m_select_B_s.set_vector(&m_B_s);
        m_A_s.load(in);
        m_rank_A_s.load(in);
        m_rank_A_s.set_vector(&m_A_s);
        m_select_A_s.load(in);
        m_select_A_s.set_vector(&m_A_s);
        m_B_o.load(in);
        m_rank_B_o.load(in);
        m_rank_B_o.set_vector(&m_B_o);
        m_select_B_o.load(in);
        m_select_B_o.set_vector(&m_B_o);
        m_A_o.load(in);
        m_rank_A_o.load(in);
        m_rank_A_o.set_vector(&m_A_o);
        m_select_A_o.load(in);
        m_select_A_o.set_vector(&m_A_o);
    }


    uint64_t size_in_bytes() {
        return sdsl::size_in_bytes(m_F_s) + sdsl::size_in_bytes(m_A_s) + sdsl::size_in_bytes(m_B_s)
               + sdsl::size_in_bytes(m_A_o) + sdsl::size_in_bytes(m_B_o);
    };


    inline uint64_t get_As(uint64_t i) {
        //cout << " get_As >> m_rank_B_s = " << i+1 << endl; fflush(stdout);
        //cout << m_rank_B_s(i+1) << endl;
        return m_select_A_s(m_rank_B_s(i + 1) + !m_B_s[i]);
    }

    inline uint64_t get_Ao(uint64_t i) {
        return m_select_A_o(m_rank_B_o(i + 1) + !m_B_o[i]);
    }

    inline uint64_t rank_Co(uint64_t i, uint64_t c) {
        return m_F_s.range_search_2d(1, i, get_As(c) + 1, get_As(c + 1), false).first;
    }

    inline uint64_t rank_Cs(uint64_t i, uint64_t c) {
        return m_F_s.range_search_2d(get_Ao(c) + 1, get_Ao(c + 1), 1, i, false).first;
    }

    inline uint64_t range_next_value_Co(uint64_t s, uint64_t e, uint64_t c) {
        //cout << " > m_F_s.rel_min_obj_maj(" << s << "," << e << "," << get_Ao(c)+1 << ")" << endl;
        //char d;
        //cin >> d;
        uint64_t t = m_F_s.rel_min_obj_maj(s, e, get_Ao(c) + 1);
        //cout << "get_Ao(" << c << ")=" << get_Ao(c) << endl;
        //cout << " > range_next_value_Co t = " << t << " m_n_p=" << m_n_p << endl;
        //cout << " > m_F_s.rel_min_obj_maj(" << s << "," << e << "," << get_Ao(c)+1 << ") = " << t << endl;

        if (t > m_n_p) return 0;
        //cout << "El rank final hasta " << t << " en Co da = " << m_rank_A_o(t) << endl;
        return m_select_B_o(m_rank_A_o(t));
    }

    inline uint64_t range_next_value_Cs(uint64_t s, uint64_t e, uint64_t c) {
        //uint64_t t = m_F_s.range_next_value(get_Ao(c)+1, s, e);
        uint64_t t = m_F_s.range_next_value(get_As(c) + 1, s, e);
        if (!t) return 0;
        //cout << "t=" << t << endl;
        //cout << "El rank final hasta " << t << " en Cs da = " << m_rank_A_s(t) << endl;
        return m_select_B_s(m_rank_A_s(t));
    }

    inline uint64_t nTriples() const {
        return m_n_p;
    }

    // Largest subject/object identifier (the universe size U used at construction)
    inline uint64_t universe() const {
        return m_B_s.size() - 2;
    }

    // Number of positions in [s, e] of m_F_s whose value is in [vs, ve]. Equivalent to count_range_search_2d of
    // adriangbrandon/sdsl-lite, written here with my_expand so that it also compiles with darroyue/sdsl-lite.
    inline uint64_t count_2d(uint64_t s, uint64_t e, uint64_t vs, uint64_t ve) const {
        typedef typename table_column_type::node_type node_type;
        struct item_type {
            node_type node;
            sdsl::range_type range;
            uint64_t ilb, irb; // interval of values of the node
        };
        uint64_t levels = m_F_s.max_level;
        uint64_t max_value = (levels >= 64) ? -1ULL : (1ULL << levels) - 1;
        if (ve > max_value) ve = max_value;
        if (s > e || vs > ve) return 0;
        uint64_t count = 0, rnk;
        sdsl::range_type r_left, r_right;
        std::vector<item_type> stack;
        stack.push_back(item_type{m_F_s.root(), {s, e}, 0, max_value});
        while (!stack.empty()) {
            item_type x = stack.back();
            stack.pop_back();
            if (m_F_s.is_leaf(x.node) || (vs <= x.ilb && x.irb <= ve)) {
                count += x.range[1] + 1 - x.range[0];
            } else {
                auto children = m_F_s.my_expand(x.node, x.range, r_left, r_right, rnk);
                uint64_t mid = ((x.ilb + x.irb) >> 1) + 1;
                if (r_left[1] + 1 > r_left[0] && vs < mid) {
                    stack.push_back(item_type{children[0], r_left, x.ilb, mid - 1});
                }
                if (r_right[1] + 1 > r_right[0] && ve >= mid) {
                    stack.push_back(item_type{children[1], r_right, mid, x.irb});
                }
            }
        }
        return count;
    }

    // Number of points (triples) in positions [s, e] of C_s and values [vs, ve] of m_F_s
    inline uint64_t count_points(uint64_t s, uint64_t e, uint64_t vs, uint64_t ve) const {
        if (s > e || vs > ve) return 0;
        return count_2d(s, e, vs, ve);
    }

    /**********************************/
    // Functions for the VEO estimators. Positions of C_s follow the object order and the values of m_F_s are
    // positions in the subject order, so an interval of C_s is a range of positions and an interval of C_o is a
    // range of values.

    // Number of triples in positions [s, e] of C_s whose subject is smaller than v
    inline uint64_t count_S_lt(uint64_t s, uint64_t e, uint64_t v) {
        if (s > e) return 0;
        uint64_t t = (v > universe()) ? m_n_p : get_As(v); // subjects < v are the values [1, t]
        if (t == 0) return 0;
        if (s == 1 && e == m_n_p) return t;
        return count_2d(s, e, 1, t);
    }

    // Number of triples in values [s, e] of C_o whose object is smaller than v
    inline uint64_t count_O_lt(uint64_t s, uint64_t e, uint64_t v) {
        if (s > e) return 0;
        uint64_t t = (v > universe()) ? m_n_p : get_Ao(v); // objects < v are the positions [1, t]
        if (t == 0) return 0;
        if (s == 1 && e == m_n_p) return t;
        return count_2d(1, t, s, e);
    }

    // Subject of the q-th triple (0-based, in subject order) among positions [s, e] of C_s
    inline uint64_t kth_S(uint64_t s, uint64_t e, uint64_t q) {
        uint64_t t;
        if (s == 1 && e == m_n_p) {
            t = q + 1;
        } else {
            auto v = m_F_s.root();
            sdsl::range_type r = {s, e}, r_left, r_right;
            uint64_t rnk;
            while (!m_F_s.is_leaf(v)) {
                auto children = m_F_s.my_expand(v, r, r_left, r_right, rnk);
                uint64_t n_left = r_left[1] + 1 - r_left[0];
                if (q < n_left) {
                    v = children[0];
                    r = r_left;
                } else {
                    q -= n_left;
                    v = children[1];
                    r = r_right;
                }
            }
            t = m_F_s.sym(v);
        }
        return m_select_B_s(m_rank_A_s(t));
    }

    // Object of the q-th triple (0-based, in object order) among values [s, e] of C_o
    inline uint64_t kth_O(uint64_t s, uint64_t e, uint64_t q) {
        uint64_t t;
        if (s == 1 && e == m_n_p) {
            t = q + 1;
        } else {
            // smallest position t such that [1, t] contains q+1 points with values in [s, e]
            uint64_t lo = 1, hi = m_n_p;
            while (lo < hi) {
                uint64_t mid = (lo + hi) / 2;
                if (count_2d(1, mid, s, e) > q) hi = mid;
                else lo = mid + 1;
            }
            t = lo;
        }
        return m_select_B_o(m_rank_A_o(t));
    }

};

#endif

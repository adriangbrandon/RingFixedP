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
//#include "bwt.hpp"
//#include "bwt_interval.hpp"
#include <configuration.hpp>

//#include <stdio.h>
//#include <stdlib.h>
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
    sub_ring(std::vector <std::pair<uint64_t, uint64_t>> &D, uint64_t U) {
        std::vector <uint64_t> C_s, C_o;
        std::vector <uint64_t> _A_s(U + 2, 0), _A_o(U + 2, 0);

        sdsl::int_vector<> _F_s(D.size() + 2);

        cout << D.size() << "  ";
        fflush(stdout);
        m_n_p = D.size();
        // First sort using the first component of the pairs
        sort(D.begin(), D.end(), [](const std::pair<uint64_t, uint64_t> &a,
                const std::pair<uint64_t, uint64_t> &b) { return a.first < b.first; });
        // Now build C_o.
        C_o.push_back(0); //dummy value
        for (uint64_t i = 0; i < D.size(); i++)
            C_o.push_back(D[i].second);

        // Then, stable sort using the second component of the pairs
        stable_sort(D.begin(), D.end(), [](const std::pair<uint64_t, uint64_t> &a,
                const std::pair<uint64_t, uint64_t> &b) { return a.second < b.second; });
        // Now build C_s
        C_s.push_back(0); //dummy value
        for (uint64_t i = 0; i < D.size(); i++)
            C_s.push_back(D[i].first);

        // Now, compute arrays m_A_s and m_A_o
        sdsl::bit_vector _B_s(U + 2, 0), _B_o(U + 2, 0);
        for (uint64_t i = 0; i < D.size(); ++i) {
            _A_s[D[i].first]++;
            _B_s[D[i].first] = 1;
            _A_o[D[i].second]++;
            _B_o[D[i].second] = 1;
        }
        _B_s[U + 1] = 1;
        _B_o[U + 1] = 1;

        uint64_t t = 0, t_p;
//            cout << "m_A_s" << endl;
        for (uint64_t i = 1; i < U + 1; ++i) {
            t_p = _A_s[i];
            _A_s[i] = t;
//                cout << _A_s[i] << endl;
            t += t_p;
        }

        t = 0;
        for (uint64_t i = 1; i < U + 1; ++i) {
            t_p = _A_o[i];
            _A_o[i] = t;
            t += t_p;
        }

        // Luego, calcular m_F_s usando LF
        std::vector <uint64_t> v_aux(U + 2, 0);
        _F_s[0] = 0; // dummy value
        //cout << "m_F_s" << endl;
        for (uint64_t i = 1; i < C_s.size(); i++) {
            v_aux[C_s[i]]++;
            _F_s[i] = _A_s[C_s[i]] + v_aux[C_s[i]];
            //cout << " > m_F_s[" << i << "] = " << _m_F_s[i] << endl;
        }
        //cout << "*******" << endl;

        // Luego construir la WM para m_F_s
        util::bit_compress(_F_s);
        construct_im(m_F_s, _F_s);

//            cout << "U="<<U <<endl;
//            cout << "n="<<D.size() << endl;
        {
            sdsl::bit_vector a(D.size() + 1, 0);

            for (uint64_t i = 0; i < _A_s.size(); i++) {
//                cout << i << " " << i+_A_s[i] << endl;
                if (_B_s[i])
                    a[_A_s[i]] = 1;
            }
            a[D.size()] = 1;
/*
                cout << "As" << endl;
                for (uint64_t i = 0; i < a.size(); i++) {
                    cout << a[i];
                }
                cout << endl; 
                cout << "Bs" << endl;
                for (uint64_t i = 0; i < m_B_s.size(); i++) {
                    cout << m_B_s[i];
                }
                cout << endl;
*/
            m_A_s = sd_vector<>(a);
            m_B_s = sd_vector<>(_B_s);
            util::init_support(m_select_A_s, &m_A_s);
            util::init_support(m_rank_A_s, &m_A_s);
            util::init_support(m_select_B_s, &m_B_s);
            util::init_support(m_rank_B_s, &m_B_s);
        }

        {
            sdsl::bit_vector a(D.size() + 1, 0);

            for (uint64_t i = 0; i < _A_o.size(); i++) {
//                cout << i << " " << i+_A_s[i] << endl;
                if (_B_o[i])
                    a[_A_o[i]] = 1;
            }
            a[D.size()] = 1;
/*
                cout << "Ao" << endl;
                for (uint64_t i = 0; i < a.size(); i++) {
                    cout << a[i];
                }
                cout << endl;
                cout << "Bo" << endl;
                for (uint64_t i = 0; i < m_B_o.size(); i++) {
                    cout << m_B_o[i];
                }
                cout << endl;
*/
            m_A_o = sd_vector<>(a);
            m_B_o = sd_vector<>(_B_o);
            util::init_support(m_rank_A_o, &m_A_o);
            util::init_support(m_select_A_o, &m_A_o);
            util::init_support(m_rank_B_o, &m_B_o);
            util::init_support(m_select_B_o, &m_B_o);
        }

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

};

#endif

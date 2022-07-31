/***
BSD 2-Clause License

Copyright (c) 2018, Adrián
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/


//
// Created by Adrián on 21/7/22.
//

#ifndef RING_LTJ_ITERATOR_HPP
#define RING_LTJ_ITERATOR_HPP

#define VERBOSE 0

namespace ring {

    template<class ring_t, class var_t, class cons_t>
    class ltj_iterator {

    public:
        typedef cons_t value_type;
        typedef var_t var_type;
        typedef ring_t ring_type;
        typedef uint64_t size_type;
        //enum state_type {s, p, o};
        //std::vector<value_type> leap_result_type;

    private:
        const triple_pattern *m_ptr_triple_pattern;
        ring_type *m_ptr_ring; //TODO: should be const
        bwt_interval m_i_s;
        bwt_interval m_i_p;
        bwt_interval m_i_o;
        value_type m_cur_s;
        value_type m_cur_p;
        value_type m_cur_o;
        bool m_is_empty = false;
        //std::stack<state_type> m_states;


        void copy(const ltj_iterator &o) {
            m_ptr_triple_pattern = o.m_ptr_triple_pattern;
            m_ptr_ring = o.m_ptr_ring;
            m_i_s = o.m_i_s;
            m_i_p = o.m_i_p;
            m_i_o = o.m_i_o;
            m_cur_s = o.m_cur_s;
            m_cur_p = o.m_cur_p;
            m_cur_o = o.m_cur_o;
            m_is_empty = o.m_is_empty;
        }

        inline bool is_variable_subject(var_type var) {
            return m_ptr_triple_pattern->term_s.is_variable && var == m_ptr_triple_pattern->term_s.value;
        }

        inline bool is_variable_predicate(var_type var) {
            return m_ptr_triple_pattern->term_p.is_variable && var == m_ptr_triple_pattern->term_p.value;
        }

        inline bool is_variable_object(var_type var) {
            return m_ptr_triple_pattern->term_o.is_variable && var == m_ptr_triple_pattern->term_o.value;
        }

    public:
        const bool &is_empty = m_is_empty;
        const bwt_interval &i_s = m_i_s;
        const bwt_interval &i_p = m_i_p;
        const bwt_interval &i_o = m_i_o;
        const value_type &cur_s = m_cur_s;
        const value_type &cur_p = m_cur_p;
        const value_type &cur_o = m_cur_o;

        ltj_iterator() = default;

        ltj_iterator(const triple_pattern *triple, ring_type *ring) {
            m_ptr_triple_pattern = triple;
            m_ptr_ring = ring;
            m_cur_s = -1;
            m_cur_p = -1;
            m_cur_o = -1;
            m_i_p = m_ptr_ring->open_POS();
            //Interval in P
            if (m_ptr_triple_pattern->term_p.value >= m_ptr_ring->nP()) {
                m_is_empty = true;
                return;
            }
            m_cur_p = m_ptr_triple_pattern->term_p.value;
            m_i_s = m_i_o = m_ptr_ring->down_P(m_ptr_triple_pattern->term_p.value);

        }

        //! Copy constructor
        ltj_iterator(const ltj_iterator &o) {
            copy(o);
        }

        //! Move constructor
        ltj_iterator(ltj_iterator &&o) {
            *this = std::move(o);
        }

        //! Copy Operator=
        ltj_iterator &operator=(const ltj_iterator &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        //! Move Operator=
        ltj_iterator &operator=(ltj_iterator &&o) {
            if (this != &o) {
                m_ptr_triple_pattern = std::move(o.m_ptr_triple_pattern);
                m_ptr_ring = std::move(o.m_ptr_ring);
                m_i_s = std::move(o.m_i_s);
                m_i_p = std::move(o.m_i_p);
                m_i_o = std::move(o.m_i_o);
                m_cur_s = o.m_cur_s;
                m_cur_p = o.m_cur_p;
                m_cur_o = o.m_cur_o;
                m_is_empty = o.m_is_empty;
            }
            return *this;
        }

        void swap(ltj_iterator &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            std::swap(m_ptr_triple_pattern, o.m_ptr_triple_pattern);
            std::swap(m_ptr_ring, o.m_ptr_ring);
            m_i_s.swap(o.m_i_s);
            m_i_o.swap(o.m_i_o);
            m_i_p.swap(o.m_i_p);
            std::swap(m_cur_s, o.m_cur_s);
            std::swap(m_cur_p, o.m_cur_p);
            std::swap(m_cur_o, o.m_cur_o);
            std::swap(m_is_empty, o.m_is_empty);
        }

        void down(var_type var, size_type c) { //Go down in the trie
            if (is_variable_subject(var)) {
                if (m_cur_o != -1 && m_cur_p != -1){
#if VERBOSE
                    std::cout << "Nothing to do" << std::endl;
#endif
                    return;
                }
                m_i_o = m_ptr_ring->down_P_S(m_i_s, c);
                m_cur_s = c;
            } else if (is_variable_object(var)) {
                if (m_cur_s != -1 && m_cur_p != -1){
#if VERBOSE
                    std::cout << "Nothing to do" << std::endl;
#endif
                    return;
                }
#if VERBOSE
                std::cout << "down_P_O" << std::endl;
#endif
                m_i_s = m_ptr_ring->down_P_O(m_i_o, m_cur_p, c);
                m_cur_o = c;
            }

        };


        void up(var_type var) { //Go up in the trie
            if (is_variable_subject(var)) {
                m_cur_s = -1;
#if VERBOSE
                std::cout << "Up in S" << std::endl;
#endif
            }else if (is_variable_object(var)) {
                m_cur_o = -1;
#if VERBOSE
                std::cout << "Up in O" << std::endl;
#endif
            }

        };

        value_type leap(var_type var) { //Return the minimum in the range
            //0. Which term of our triple pattern is var
            if (is_variable_subject(var)) {
                //1. We have to go down through s
                if (m_cur_p != -1 && m_cur_o != -1) {
                    //PO->S
#if VERBOSE
                    std::cout << "min_S_in_PO" << std::endl;
#endif
                    return m_ptr_ring->min_S_in_PO(m_i_s);
                } else if (m_cur_p != -1) {
                    //P->S
#if VERBOSE
                    std::cout << "min_S_in_P" << std::endl;
#endif
                    return m_ptr_ring->min_S_in_P(m_i_s);
                }
            } else if (is_variable_object(var)) {
                //1. We have to go down in the trie of o
                if (m_cur_s != -1 && m_cur_p != -1) {
                    //SP->O
#if VERBOSE
                    std::cout << "min_O_in_SP" << std::endl;
#endif
                    return m_ptr_ring->min_O_in_PS(m_i_o);
                } else if (m_cur_p != -1) {
                    //P->O
#if VERBOSE
                    std::cout << "min_O_in_P" << std::endl;
#endif
                    return m_ptr_ring->min_O_in_P(m_i_o);
                }
            }
            return 0;
        };

        value_type leap(var_type var, size_type c) { //Return the next value greater or equal than c in the range
            //0. Which term of our triple pattern is var
            if (is_variable_subject(var)) {
                //1. We have to go down through s
                if (m_cur_p != -1 && m_cur_o != -1) {
                    //PO->S
#if VERBOSE
                    std::cout << "next_S_in_PO" << std::endl;
#endif
                    return m_ptr_ring->next_S_in_PO(m_i_s, c);
                } else if (m_cur_p != -1) {
                    //P->S
#if VERBOSE
                    std::cout << "next_S_in_P" << std::endl;
#endif
                    return m_ptr_ring->next_S_in_P(m_i_s, c);
                }
            } else if (is_variable_object(var)) {
                //1. We have to go down in the trie of o
                if (m_cur_s != -1 && m_cur_p != -1) {
                    //SP->O
#if VERBOSE
                    std::cout << "next_O_in_SP" << std::endl;
#endif
                    return m_ptr_ring->next_O_in_PS(m_i_o, c);
                } else if (m_cur_p != -1) {
                    //P->O
#if VERBOSE
                    std::cout << "next_O_in_P" << std::endl;
#endif
                    return m_ptr_ring->next_O_in_P(m_i_o, c);
                }
            }
            return 0;
        }

        bool in_last_level(){
            return (m_cur_o !=-1 && m_cur_p != -1) || (m_cur_s !=-1 && m_cur_p != -1)
                    || (m_cur_o !=-1 && m_cur_s != -1);
        }

        //Solo funciona en último nivel, en otro caso habría que reajustar
        std::vector<uint64_t> seek_all(var_type var){
            if (is_variable_subject(var)){
                return m_ptr_ring->all_S_in_range(m_i_s);
            }else if (is_variable_object(var)){
                return m_ptr_ring->all_O_in_range(m_i_o);
            }
            return std::vector<uint64_t>();
        }
    };

}

#endif //RING_LTJ_ITERATOR_HPP

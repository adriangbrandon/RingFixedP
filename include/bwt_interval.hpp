/*
 * bwt_interval.hpp
 * Copyright (C) 2020 Author removed for double-blind evaluation
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

#ifndef BWT_INTERVAL
#define BWT_INTERVAL



namespace ring {
    class bwt_interval {
        uint64_t l;
        uint64_t r;
        uint64_t cur_val;  // current value within the interval
        uint64_t cur_rank;
        uint64_t p_value;

    private:

        void copy(const bwt_interval &o){
            l = o.l;
            r = o.r;
            cur_val = o.cur_val;
            cur_rank = o.cur_rank;
            p_value = o.p_value;
        }

    public:
        bwt_interval() = default;


        bwt_interval(uint64_t _l, uint64_t _r)
        {
            l = _l;
            r = _r;
            cur_val = (uint64_t)-1;
            cur_rank = (uint64_t)-1;
            p_value = (uint64_t)-1;
        }

        bwt_interval(uint64_t _l, uint64_t _r, uint64_t _p_value)
        {
            l = _l;
            r = _r;
            cur_val = (uint64_t)-1;
            cur_rank = (uint64_t)-1;
            p_value = _p_value;
        }


        //! Copy constructor
        bwt_interval(const bwt_interval &o) {
            copy(o);
        }

        //! Move constructor
        bwt_interval(bwt_interval &&o) {
            *this = std::move(o);
        }

        //! Copy Operator=
        bwt_interval &operator=(const bwt_interval &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }

        //! Move Operator=
        bwt_interval &operator=(bwt_interval &&o) {
            if (this != &o) {
                l = o.l;
                r = o.r;
                cur_val = o.cur_val;
                cur_rank = o.cur_rank;
                p_value = o.p_value;
            }
            return *this;
        }

        void swap(bwt_interval &o) {
            std::swap(l, o.l);
            std::swap(r, o.r);
            std::swap(cur_val, o.cur_val);
            std::swap(cur_rank, o.cur_rank);
            std::swap(p_value, o.p_value);
        }


        inline uint64_t get_cur_value() {
            return cur_val;
        }

        inline uint64_t get_p_value() const
        {
            return p_value;
        }

        inline uint64_t left() {
            return l;
        }

        inline uint64_t right() {
            return r;
        }

        inline uint64_t size() const {
            return r - l + 1;
        }
    };
}

#endif

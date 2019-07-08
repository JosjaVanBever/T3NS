/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018-2019 Josja Van Bever <Josja.VanBever@UGent.be>
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "cpp_interface.h"
#include "symsecs.h"


template <class T>
class Pair {
    public:
        // constructor
        Pair(T val1, T val2) : data(val1,val2);
        // access operators
        T & operator[](int ind) {return data[ind]; };       
        T const & operator[](int ind) const {return data[ind]; };
    private:
        T data[2];
};

class SymsecMatcher {
    public:
        // constructor and destructor
        SymsecMatcher(int max_size=0) : size(0), max_size(max_size) {
            result = (Pair<int> *) safe_malloc(max_size, Pair<int>); }
//            result = (int (*)[2]) safe_malloc(2 * max_size, int); }
        ~ SymsecMatcher() { free(result); }

        // get the current results
        int get_size() const { return size; };
        const Pair<int> * get_result() const { return result; };
//        const int * get_result()[2] const { return result; };

        void set_size(int new_size) {
            if (new_size > max_size) {
                // WARNING: MEMORY LEAK MIGHT BE POSSIBLE?
                this->result = (Pair<int> *) realloc(result, new_size);
            }
            size = new_size;
        }

        // @TEST
        void set_test_result(int startval, int size) {
            this->size = size;
            for (int i=0; i<size; i++)  {
                result[i][0]  = (startval != 0)? startval+i : 0;
                result[i][1]  = (startval != 0)? -startval-i : 0;
            }
        }

        // // calculate the new results
        // // The results are guaranteed to be sorted in the irrep index
        // // of the first argument.
        // void set_matching_symsec_indices(const TensorInfo * a,
        //     const TensorInfo * b, int leg) const;
        // void set_matching_symsec_indices(const struct symsecs * a,
        //     const struct symsecs * b) const;
    private:
        // data fields
        int size;  // effective space currently used
        int max_size;  // allocated space
        Pair<int> * result;
//        int (*result)[2];  // pointer to array of int pairs
};
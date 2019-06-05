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

#ifndef OVERLAP_CALCULATOR_H
#define OVERLAP_CALCULATOR_H

#include "cpp_interface.h"
#include "bookkeeper.h"
#include "tensor_info.h"
#include <iostream>

// The overlap calculator does not allocate any of the T3NS data operated on.
// It does allocate its own temporary results.
class OverlapCalculator {
    public:
        // constructor and destructor
        OverlapCalculator(T3NSfill* opt, T3NSfill* ref);
        ~OverlapCalculator();

        // @TEST
        int get_result();
    private:
        // Main data:
        struct TensorInfoPair * tensorpairs;
        int nr_tensorpairs;  // length of tensorpairs

        // Help data:
        struct bookkeeper * opt_bookie;
        struct bookkeeper * ref_bookie;

        // @TEST
        int test;
};

#endif

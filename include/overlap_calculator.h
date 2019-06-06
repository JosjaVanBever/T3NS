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

#include "overlap_object.h"

// interface
#include "cpp_interface.h"
// C libraries
#include "bookkeeper.h"
// C++ libraries
#include "siteTensor.h"
#include "tensor_info.h"

// external libraries
#include <iostream>
#include <stdlib.h>
#include <time.h>

#ifndef OVERLAP_OBJECT_H
#error 'OverlapObject not defined!'
#endif

// The overlap calculator does not allocate any of the T3NS data operated on.
// It does allocate its own temporary results.
class OverlapCalculator {
    public:
        // constructor and destructor
        OverlapCalculator(const T3NSfill * opt, const T3NSfill * ref,
                const struct network * netw);
        ~OverlapCalculator();

        // @TEST: temporary public
        // -> get the OO link attached to who and pointing to other
        int get_internal_link(int who, int other, OverlapObjectLink * result);
        // -> get the OO links attached to who but not pointing to exclude
        int get_external_links(int who, int exclude, OverlapObjectLink * result);

        // @TEST
        int get_result();
    private:
        // Main data:
        struct TensorInfoPair * tensorpairs;
            int nr_tensorpairs;  // length of tensorpairs
        OverlapObject * overlaps;
            int nr_overlaps;  // length of overlaps
        // remark: there are also overlaps living on a bond to -1!

        // Help data:
        const struct bookkeeper * opt_bookie;
        const struct bookkeeper * ref_bookie;
        const struct network * network;

        // Help functions:
        // Get the appropriate OO links from 'who' with respect to 'other',
        //   using the criterium 'check' and performing the action
        //   question(to, who) to the matches 'to' with the criterium.
        int get_links(int who, int other, bool (*check)(int,int,int),
                void (*question)(int,int), OverlapObjectLink * result);
};

// External help functions:
// -> is other in (from,to)?
bool found_other(int other, int from, int to);
// -> is other not in (from, to)?
bool avoid_other(int other, int from, int to);

#endif

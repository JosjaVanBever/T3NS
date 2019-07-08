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


#include "symsec_matcher.h"
#include "symsecs.h"
#include "tensor_info.h"


void SymsecMatcher::set_size(int new_size)
{
    // change max_size if necessary
    if (new_size > max_size) {
        // WARNING: MEMORY LEAK MIGHT BE POSSIBLE?
        this->result = (Pair<int> *) realloc(result, new_size); }
    // change size
    size = new_size;
}


void SymsecMatcher::set_matching_symsec_indices(const TensorInfo * A,
            const TensorInfo * B, int leg)
{
    set_matching_symsec_indices(A->get_sym(leg), B->get_sym(leg));
}

void SymsecMatcher::set_matching_symsec_indices(const struct symsecs * a,
            const struct symsecs * b)
{
    int index = 0;
    for (int i = 0; i < a->nrSecs; i++) {
        int j = search_symsec((a->irreps)[i], b);
        if (j != -1) {
            (this->result)[index] = {i,j};
            index++;
        }
    }
    size = index;
}

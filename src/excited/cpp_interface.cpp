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

#ifndef CPP_INTERFACE
#warning 'cpp_interface.h not included'
#endif

// void init_overlap_calculator(int test, OverlapCalculator ** result) {
// 	*result = new OverlapCalculator(test);
// }

// #ifndef TESTT3NSFILL
// #warning 'T3NSfill was not defined'
// #endif

void init_overlap_calculator(T3NSfill* opt_t3ns, OverlapCalculator ** result) {
    *result = new OverlapCalculator(opt_t3ns);
}

int get_result(OverlapCalculator* calc) {
	return calc->get_result();
}

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


#ifndef TWOSITEOVERLAPCALCULATOR_H
#define TWOSITEOVERLAPCALCULATOR_H


#include "overlap_calculator.h"
class TensorInfo;
class OverlapObject;

class TwoSiteOverlapCalculator : public OverlapCalculator {
	public:
		// constructor and destructor
        TwoSiteOverlapCalculator(const T3NSfill * opt, const T3NSfill * ref,
                const struct network * netw);
        ~TwoSiteOverlapCalculator();

        // WE WILL ALLOCATE; WARNING FOR COPYING!
	
		// @TEST
        int perform_testing();
        int perform_set_OO_to_contraction_test();
        // void can_you_find_me() {printf("Here you are!");};

	// private:
		void set_OO_to_contraction(const TensorInfo * A, const TensorInfo * B,
				int leg, OverlapObject * OO);
		// storage of intermediate results
		TensorInfo tensmem[2]; // tensor memory
		// backup tensor memory for branching tensor calculations
		TensorInfo tensbmem[2];
		// indices of the last optimization center used
		// the indices are stored such that last_optimized ...
		int last_optimized[2];
};


#endif

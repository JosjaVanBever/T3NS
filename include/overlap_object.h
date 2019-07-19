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

#ifndef OVERLAP_OBJECT_H
#define OVERLAP_OBJECT_H

#include "cpp_interface.h"
#include "sparseblocks.h"
#include "symsecs.h"
#include <stdio.h>

class SymsecMatcher;

// IDEA: DO NOT REMEMBER REF AND OPT BUT TAKE IT AS INPUT

// Usage: ones created, first renew the block layout with a solved
// symsecMatcher before performing block calculations that take block info.
class OverlapObject {
	public:
		// constructors and destructor
		// @param:
		//   ref = symmetry sectors of the bond in the reference T3NS
		//   opt = symmetry sectors of the bond in the optimizing T3NS
		//   [opt_dim] = current/bond dimension of opt
		// @remark: if opt_dim is set to the bond dimension, no reallocation
		//   will be needed. Remark that 2000**2 * 100 * 8 bytes is about 3Gb.
		OverlapObject(struct symsecs * ref, struct symsecs * opt,
				int opt_dim = -1);  // <0 => opt->totaldims
		OverlapObject(const OverlapObject &);  // Copy constructor
		~OverlapObject();  // destructor

		// Assignment operator
		OverlapObject& operator=(const OverlapObject &);
		
		// renew the beginblocks and tel array
		// @param: match should contain the matches between ref and opt;
		//   no assumptations are made about the ordering of matching symsecs
		void renew_block_layout(const SymsecMatcher * match, bool set_zero);

		// reallocate the element array explicitely
		// @param: opt_dim = total/bond dimension for the optimizing bond
		void reallocate_optdim(int opt_dim);

		// get or set the reference or optimizing symsec
		struct symsecs * get_ref() { return ref; }
		struct symsecs * get_opt() { return opt; }
		//void set_ref(struct symsecs * new_ref) { ref = new_ref; } REALLOCATE!
		//void set_opt(struct symsecs * new_opt) { opt = new_opt; } DO NOT DO THIS EITHER!
		// get the dimensions of block i (indexing based on ref)
		int get_ldim(int i) const { return ldim[i]; }  // leading, based on ref
		int get_sdim(int i) const { return sdim[i]; }  // second, based on opt
		int get_nr_blocks() const { return ref->nrSecs; }  // number of blocks
		// get the help variables
		int get_usedsize() const { return usedsize; }
		int get_allocsize() const { return allocsize; }

	private:
		// // help function for reallocation of the element array
		// void reallocate_elements(int elements);
		// help function for copying
		void copy (const OverlapObject & copy);

		// Main data:
		// Symmetry sectors in the reference and optimizing bond
		//   to wich the OO is associated
		struct symsecs * ref, * opt;
		// Actual elements of the OverlapObject.
		//   The block corresponding to the i'th irrep of ref starts at   ==> label with opt instead
		//   blocks.tel[blocks.beginblock[i]], has leading dimension
		//   ldim[i] and second dimension sdim[i].
		//   If blocks.beginblock[i] == -1, the i'th block is not present.  ==> can be used as check
		struct sparseblocks blocks;  // begin blocks and element array
		int * ldim;     // leading dimension array -> from reference
		int * sdim;     // second dimension array -> from optimizing

		// Help data:
		// size of tel array used to effectively store elements
		int usedsize;
		// allocated size of the tel array; should be >= usedsize
		int allocsize;
};


void print_overlap_object(const struct bookkeeper * ref_bookie,
		const struct bookkeeper * opt_bookie, OverlapObject * OO);


typedef struct OverlapObjectLink {
	OverlapObject* OO;
	int leg;  // index of the leg to which the OOlink is associated
} OverlapObjectLink;


#endif

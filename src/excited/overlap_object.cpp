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


#include "overlap_object.h"
#include "symsec_matcher.h"
#include "sparseblocks.h"
#include "cpp_tools.h"

OverlapObject::OverlapObject(struct symsecs * ref, struct symsecs * opt,
			int opt_dim) : ref(ref), opt(opt)
{
	// redefine opt_dim if smaller than 0
	if (opt_dim < 0) { opt_dim = opt->totaldims; }

	// help variable
	int nrBlocks = ref->nrSecs;

	// initialize the sparseblocks
	init_memory_sparseblocks(&blocks, nrBlocks, ref->totaldims * opt_dim);
	allocsize = ref->totaldims * opt_dim;
	usedsize = 0;

	// allocate the ldim and sdim array
	ldim = (int *) safe_malloc(nrBlocks, int);
	sdim = (int *) safe_malloc(nrBlocks, int);
}


// Copy constructor
OverlapObject::OverlapObject(const OverlapObject & copy) : ref(copy.ref),
		opt(copy.opt), usedsize(copy.usedsize), allocsize(copy.allocsize)
{
	this->copy(copy);
}

// Destructor
OverlapObject::~OverlapObject()
{
	// destroy_sparseblocks(struct sparseblocks * blocks);
	destroy_sparseblocks(&blocks);
	safe_free(ldim); safe_free(sdim);
}

// Assignment operator
OverlapObject& OverlapObject::operator=(const OverlapObject & rhs)
{
	ref = rhs.ref; opt = rhs.opt;
	usedsize = rhs.usedsize; allocsize = rhs.allocsize;
	this->copy(rhs);
	return *this;
}

// help function to copy
void OverlapObject::copy (const OverlapObject & copy)
{
	// help variable
	int nr_blocks = copy.get_nr_blocks();
	// sparseblocks
	deep_copy_sparseblocks(&(this->blocks), &(copy.blocks), nr_blocks);
	// ldim and sdim
	ldim = (int *) safe_malloc(nr_blocks, int);
	sdim = (int *) safe_malloc(nr_blocks, int);
	for (int i=0; i<nr_blocks; i++) {
		ldim[i] = copy.ldim[i];
		sdim[i] = copy.sdim[i];
	}
}

// void OverlapObject::set_ref(struct symsecs * new_ref)
// {
// 	// help variable
// 	int nrBlocks = new_ref->nrSecs;
// 	// reallocate if necessary
// 	if (nrBlocks > get_nr_blocks()) {
// 		reallocate<int>(ldim, nrBlocks);
// 		reallocate<int>(sdim, nrBlocks);
// 		reallocate<int>(blocks.beginblock, nrBlocks);
// 	}
// 	// values and length of tel array are not meaningfull anymore
// 	blocks.beginblock[nrBlocks] = 0;

// 	// set the ref variable
// 	this->ref = new_ref;
// }

// renew the beginblocks and tel array
// match should contain the matches between ref and opt
// no assumptations are made about the ordering of matching symsecs
void OverlapObject::renew_block_layout(const SymsecMatcher * match,
		bool set_zero)
{
	// help variables
	usedsize = 0;
	int nrBlocks = ref->nrSecs;

	// reset the beginblock and dimension arrays
	for (int i=0; i<nrBlocks; i++) {
		blocks.beginblock[i] = -1;  // default value
	 	ldim[i] = sdim[i] = 0; }    // start empty
	blocks.beginblock[nrBlocks] = -1;

	// fill the beginblock and dimension arrays
	for (int j=0; j<match->get_size(); j++) {
		// get the relevant indices
	 	int ref_index = match->get_result()[j][0];
	 	int opt_index = match->get_result()[j][1];
	 	int index = opt_index;  // index used for labeling blocks   @CHANGED

	 	assert(index < nrBlocks);

		// set the begin and dims of the block
	 	ldim[index] = ref->dims[ref_index];
		sdim[index] = opt->dims[opt_index];
		// if the block is non empty
		if (ldim[index] && sdim[index] ) {   // true if != 0   @CHANGED
			blocks.beginblock[index] = usedsize;
			// increase the usedsize
	 		usedsize += ldim[index] * sdim[index];
		}	
	}
	blocks.beginblock[nrBlocks] = usedsize;

	// reallocate tel if necessary
	if (usedsize > allocsize) {
		allocsize = 2 * usedsize;
		reallocate<EL_TYPE>(&blocks.tel, allocsize); }
	if (set_zero) {
		// fill tel up with zeros
		for (int i=0; i<usedsize; i++) {
	 		blocks.tel[i] = 0;
		}
	}
}


void OverlapObject::reallocate_optdim(int opt_dim)
{
	// reallocate_elements(ref->totaldims * opt_dim);
	reallocate<EL_TYPE>(&blocks.tel, ref->totaldims * opt_dim);
}


// void OverlapObject::reallocate_elements(int elements)
// {
// 	allocsize = elements;
// 	blocks.tel = (EL_TYPE *) realloc(blocks.tel, elements * sizeof(EL_TYPE));
// }


void print_overlap_object(const struct bookkeeper * ref_bookie,
		const struct bookkeeper * opt_bookie, OverlapObject * OO)
{
	fprintf(stdout, "Symmetries:\n");
	fprintf(stdout, "ref: ");
	print_symsecs(ref_bookie, OO->get_ref(), 0);
	fprintf(stdout, "opt: ");
	print_symsecs(opt_bookie, OO->get_opt(), 0);

	fprintf(stdout, "Blocks:\n");
	for (int i=0; i<OO->get_nr_blocks(); i++) {
		fprintf(stdout, "(%d,%d) ", OO->get_ldim(i), OO->get_sdim(i));
	}
	fprintf(stdout, "\n");

	fprintf(stdout, "Internal data:\n");
	fprintf(stdout, "nr_blocks: %d\n", OO->get_nr_blocks());
	fprintf(stdout, "usedsize: %d\n", OO->get_usedsize());
	fprintf(stdout, "allocsize: %d\n", OO->get_allocsize());
}
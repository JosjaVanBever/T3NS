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

OverlapObject::OverlapObject(struct symsecs * ref, struct symsecs * opt,
			int opt_dim) : ref(ref), opt(opt)
{
	// redefine opt_dim if smaller than 0
	if (opt_dim < 0) { opt_dim = 2 * opt->totaldims; }

	// help variable
	int nrBlocks = ref->nrSecs;

	// initialize the sparseblocks
	init_null_sparseblocks(&blocks);
	// -> allocate the tel array
	blocks.tel = (EL_TYPE *) safe_malloc(ref->totaldims * opt_dim, EL_TYPE);
	allocsize = opt_dim;
	// -> allocate beginblock
	blocks.beginblock = (int *) safe_malloc(nrBlocks, int);

	// allocate the ldim and sdim array
	ldim = (int *) safe_malloc(nrBlocks, int);
	sdim = (int *) safe_malloc(nrBlocks, int);
}


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
	for (int i=0; i<nrBlocks+1; i++) {
		blocks.beginblock[i] = -1;  // default value
		ldim[0] = sdim[0] = 0;
	}
	// fill the beginblock and dimension arrays
	for (int j=0; j<match->get_size(); j++) {
		// get the relevant indices
		int ref_index = match->get_result()[j][0];
		int opt_index = match->get_result()[j][1];
		// set the begin and dims of the block
		blocks.beginblock[ref_index] = usedsize;
		ldim[ref_index] = ref->dims[ref_index];
		sdim[opt_index] = opt->dims[opt_index];
		// increase the usedsize
		usedsize += ldim[ref_index] * sdim[opt_index];
	}
	// reallocate tel if necessary
	if (usedsize > allocsize) {
		allocsize = 2 * usedsize;
		this->reallocate_elements(allocsize); }
	if (set_zero) {
		// fill tel up with zeros
		for (int i=0; i<usedsize; i++) {
			blocks.tel[i] = 0; }
	}
}


OverlapObject::~OverlapObject()
{
	safe_free(blocks.beginblock); safe_free(blocks.tel);
	safe_free(ldim); safe_free(sdim);
}


void OverlapObject::reallocate_optdim(int opt_dim)
{
	reallocate_elements(ref->totaldims * opt_dim);
}


void OverlapObject::reallocate_elements(int elements)
{
	allocsize = elements;
	blocks.tel = (EL_TYPE *) realloc(blocks.tel, elements * sizeof(EL_TYPE));
}


void print_overlap_object(const struct bookkeeper * ref_bookie,
		const struct bookkeeper * opt_bookie, OverlapObject * OO)
{
	fprintf(stdout, "ref: ");
	print_symsecs(ref_bookie, OO->get_ref(), 0);
	fprintf(stdout, "opt: ");
	print_symsecs(opt_bookie, OO->get_opt(), 0);
}
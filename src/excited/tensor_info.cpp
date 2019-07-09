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

#include "tensor_info.h"
#include "overlap_object.h"

TensorInfo::TensorInfo(struct siteTensor * data, struct symsecs ** syms)
	: data(data)
{
	for (int i=0; i<3; i++) { this->syms[i] = syms[i]; }
}


struct symsecs * TensorInfo::get_sym(int leg) const
{
	return syms[leg];
}

void TensorInfo::get_syms(struct symsecs ** result) const
{
	for (int i=0; i<3; i++) { result[i] = syms[i]; }
}


void TensorInfo::renew_symsec_layout(const TensorInfo * ref,
		const OverlapObjectLink * OO_link)
{
	// extract the index of the leg over which is contracted
	int contracted_leg = OO_link->leg;

	// check whether the contraction is valid
	assert(OO_link->OO->get_ref() == ref->get_sym(OO_link->leg));

	// set the relevant sectors
	this->set_syms(ref->get_syms());
	this->set_sym(contracted_leg, OO_link->OO->get_opt());
}


// void renew_block_layout(const TensorInfo * ref, bool set_zero)
// {
// 	// number of blocks of the new structure
// 	int nrblocks = ref->get_data()->nrblocks;

// 	// reallocate qnumbers and beginblock if necessary
// 	if (data->nrblocks < nrblocks) {
// 		realloc(data->qnumbers, nrblocks, sizeof(QN_TYPE));
// 		realloc(data->blocks->beginblock, nrblocks, sizeof(int));
// 		nr_allocated_blocks = nrblocks;
// 	}

// 	// set number of blocks
// 	data->nrblocks = nrblocks;

// 	// fill qnumbers
// 	for (int i=0; i<data->nrblocks; i++) {
// 		data->qnumbers[i] = ref->data->qnumbers[i];
// 	}

// 	// fill beginblock
// 	int usedsize = 0;
// 	int dims[3];  // dimensions of the current block
// 	for (int i=0; i<data->nrblocks; i++) {
// 		get_block_dimensions(i, dims);
// 		data->blocks->beginblock[i] = usedsize;
// 		usedsize += dims[0] * dims[1] * dims[2];
// 	}

// 	// reallocate tel if necessary
// 	if (usedsize > allocsize) {
// 		allocsize = 2 * usedsize;
// 		this->reallocate(allocsize); }
// 	if (set_zero) {
// 		// fill tel up with zeros
// 		for (int i=0; i<usedsize; i++) {
// 			blocks.tel[i] = 0; }
// 	}
// }


void print_tensorInfo(const struct bookkeeper * keeper,
		const TensorInfo * tensor)
{
	// print the SiteTensor structure
	fprintf(stdout,"SiteTensor data:\n");
	print_siteTensor(keeper, tensor->get_data());

	// print the Symsecs structures
	fprintf(stdout,"Symsectors array:\n");
	struct symsecs * syms[3];
	tensor->get_syms(syms);
	for (int i=0; i<3; i++) {
		print_symsecs(keeper, syms[i], 0);
	}
}


void print_TensorInfoPair(const struct bookkeeper * opt_bookie,
		const struct bookkeeper * ref_bookie,
		const struct TensorInfoPair * pair)
{
	// print the opt TensorInfo
	fprintf(stdout,"TensorInfoPair with\nOPT:\n");
	print_tensorInfo(opt_bookie, &(pair->opt));

	// print the ref TensorInfo
	fprintf(stdout,"REF:\n");
	print_tensorInfo(ref_bookie, &(pair->ref));
}

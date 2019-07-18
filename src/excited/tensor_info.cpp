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
#include "cpp_interface.h"
#include "overlap_object.h"
#include <stdio.h>
#include "siteTensor.h"
#include "cpp_tools.h"


TensorInfo::TensorInfo() : data(NULL), is_physical(false), usedsize(0), allocsize(0),
		nr_allocated_blocks(0)
{
	for (int i=0; i<3; i++) { this->syms[i] = NULL; }
}

TensorInfo::TensorInfo(struct siteTensor * data, struct symsecs ** syms, bool is_psite)
		: data(data), is_physical(is_psite), nr_allocated_blocks(data->nrblocks)
{
	usedsize = allocsize = (data->blocks.beginblock != NULL)?
			data->blocks.beginblock[data->nrblocks] : 0;
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


void TensorInfo::copy_symmetry_layout(const TensorInfo * ref)
{
	// set correct type of tensor
	this->is_physical = ref->is_physical;

	// set the sectors
	for (int i=0; i<3; i++) {
		this->syms[i] = ref->syms[i];
	}

	printf("temporary solution used for symsecs!\n");
}


// void TensorInfo::renew_symsec_layout(const TensorInfo * ref,
// 		const OverlapObjectLink * OO_link)
// {
// 	// check whether the contraction is valid
// 	assert(OO_link->OO->get_ref() == ref->get_sym(OO_link->leg));

// 	// set correct type of tensor
// 	this->is_physical = ref->is_physical;

// 	// set the sectors
// 	for (int i=0; i<3; i++) {
// 		this->syms[i] = ref->syms[i];
// 	}
// }


void TensorInfo::renew_symsec_layout(const TensorInfo * ref,
		const OverlapObjectLink * OO_link)
{
	// set correct type of tensor
	this->is_physical = ref->is_physical;

	// extract the index of the leg over which is contracted
	int contracted_leg = OO_link->leg;

	// check whether the contraction is valid
	assert(OO_link->OO->get_ref() == ref->get_sym(OO_link->leg));

	// set the relevant sectors
	for (int i=0; i<3; i++) {
		this->syms[i] = (i == contracted_leg)? OO_link->OO->get_opt() :
				ref->syms[i];
	}
}


////////////////////////////////
inline int TensorInfo::get_block_dimension(int leg, int sym_index) const
{
	return (syms[leg])->dims[sym_index];
}


inline void TensorInfo::get_block_dimensions(int blocknr, int * dims) const
{
	// get the leg indices
	int ids[3];
	get_sym_indices(blocknr, ids);
	// extract the corresponding dimensions
	for (int i=0; i<3; i++) {
		dims[i] = get_block_dimension(i, ids[i]); }
}

inline int TensorInfo::get_sym_index(int blocknr, int leg) const
{
	int ids[3];
	get_sym_indices(blocknr, ids);
	return ids[leg];
}

inline void TensorInfo::get_sym_indices(int blocknr, int * result) const
{
	indexize_address(result, data->qnumbers[blocknr], syms);
}
/////////////////////////////////////

void TensorInfo::renew_block_layout(const TensorInfo * ref,
		const struct OverlapObjectLink * OO_link, bool set_zero)
{
	// number of blocks of the new structure
	int nrblocks = ref->data->nrblocks;
	// int nrblocks = OO_link->OO->get_nr_blocks();  //WRONG!!!!!!!
	// set number of blocks
	data->nrblocks = nrblocks;

	// reallocate qnumbers and beginblock if necessary
	if (nr_allocated_blocks < nrblocks) {
		reallocate<QN_TYPE>(&(data->qnumbers), nrblocks);
		reallocate<int>(&(data->blocks.beginblock), nrblocks);
		nr_allocated_blocks = nrblocks;
	}

	// fill qnumbers
	for (int i=0; i<nrblocks; i++) {
		data->qnumbers[i] = translate_qn_address(ref->data->qnumbers[i], ref->syms, this->syms);
		//data->qnumbers[i] = ref->data->qnumbers[i]; // hier translate qnumber als symsec aangepast!!!
	}

	// fill beginblock
	usedsize = 0;
	int dims[3];  // dimensions of the current block
	int ids[3];   // symmetry indices of the current block
	for (int i=0; i<nrblocks; i++) {
		get_sym_indices(i, ids);  //impliciet qnumbers bekijken!!!
		// loop over all legs
		for (int j=0; j<3; j++) {
			// set the appropriate dimension
			dims[j] = (j == OO_link->leg)? OO_link->OO->get_sdim(ids[j]) :   //WRONG!!!
					get_block_dimension(j, ids[j]); }
		// set the element in beginblock
		data->blocks.beginblock[i] = usedsize;
		// increase the usedsize
		usedsize += dims[0] * dims[1] * dims[2];

		fprintf(stdout, "beginblock[%d]: %d\n", i, data->blocks.beginblock[i]);
		fprintf(stdout, "block dimensions %d: %d %d %d\n", i, dims[0], dims[1], dims[2]);
		fprintf(stdout, "      symindices %d: %d %d %d\n", i, ids[0], ids[1], ids[2]);

	}

	// reallocate tel if necessary
	if (usedsize > allocsize) {
		allocsize = usedsize;
		reallocate<EL_TYPE>(&(data->blocks.tel), allocsize);
	}

	if (set_zero) {
		// fill tel up with zeros
		for (int i=0; i<usedsize; i++) {
			data->blocks.tel[i] = 0; }
	}
}


// void TensorInfo::renew_block_layout(const TensorInfo * ref,
// 		const struct OverlapObjectLink * OO_link, bool set_zero)
// {
// 	// number of blocks of the new structure
// 	int nrblocks = ref->data->nrblocks;
// 	// set number of blocks
// 	data->nrblocks = nrblocks;

// 	// reallocate qnumbers and beginblock if necessary
// 	if (nr_allocated_blocks < nrblocks) {
// 		reallocate<QN_TYPE>(&(data->qnumbers), nrblocks);
// 		reallocate<int>(&(data->blocks.beginblock), nrblocks);
// 		nr_allocated_blocks = nrblocks;
// 	}

// 	// fill qnumbers
// 	for (int i=0; i<nrblocks; i++) {
// 		data->qnumbers[i] = ref->data->qnumbers[i]; // hier translate qnumber
// 	}

// 	// fill beginblock
// 	usedsize = 0;
// 	int dims[3];  // dimensions of the current block
// 	int ids[3];   // symmetry indices of the current block
// 	for (int i=0; i<nrblocks; i++) {
// 		get_sym_indices(i, ids);
// 		// loop over all legs
// 		for (int j=0; j<3; j++) {
// 			// set the appropriate dimension
// 			dims[j] = (j == OO_link->leg)? OO_link->OO->get_sdim(i) :
// 					get_block_dimension(j, ids[j]); }
// 		// increase the usedsize
// 		usedsize += dims[0] * dims[1] * dims[2];
// 	}

// 	// reallocate tel if necessary
// 	if (usedsize > allocsize) {
// 		allocsize = 2 * usedsize;
// 		reallocate<EL_TYPE>(&(data->blocks.tel), allocsize);
// 	}

// 	if (set_zero) {
// 		// fill tel up with zeros
// 		for (int i=0; i<usedsize; i++) {
// 			data->blocks.tel[i] = 0; }
// 	}
// }



// void TensorInfo::renew_block_layout(const TensorInfo * ref, bool set_zero)
// {
// 	printf("inside renew_block_layout\n");
// 	fflush(stdout);

// 	// number of blocks of the new structure
// 	int nrblocks = ref->data->nrblocks;

// 	// reallocate qnumbers and beginblock if necessary
// 	if (nr_allocated_blocks < nrblocks) {
// 		reallocate<QN_TYPE>(&(data->qnumbers), nrblocks);
// 		reallocate<int>(&(data->blocks.beginblock), nrblocks);
// 		nr_allocated_blocks = nrblocks;
// 	}

// 	printf("somewhere over the rainbow earlier");
// 	fflush(stdout);

// 	// fill qnumbers
// 	for (int i=0; i<nrblocks; i++) {
// 		data->qnumbers[i] = ref->data->qnumbers[i];  // data is NULL!!!!
// 	}

// 	printf("somewhere over the rainbow");
// 	fflush(stdout);

// 	// set number of blocks
// 	data->nrblocks = nrblocks;


// 	int testdims[3];
// 	printf("\nthis:\n");
// 	for (int i=0; i<this->data->nrblocks; i++) {
// 		this->get_block_dimensions(i, testdims);
// 		printf("beginblock[%d]: %d\n", i, this->data->blocks.beginblock[i]);
// 		printf("block %d: %d %d %d\n", i, testdims[0], testdims[1], testdims[2]);
// 	}
// 	printf("\nref:\n");
// 	for (int i=0; i<ref->data->nrblocks; i++) {
// 		ref->get_block_dimensions(i, testdims);
// 		printf("beginblock[%d]: %d\n", i, ref->data->blocks.beginblock[i]);
// 		printf("block %d: %d %d %d\n", i, testdims[0], testdims[1], testdims[2]);
// 	}


// 	// fill beginblock
// 	usedsize = 0;
// 	int dims[3];  // dimensions of the current block
// 	for (int i=0; i<nrblocks; i++) {
// 		get_block_dimensions(i, dims);
// 		data->blocks.beginblock[i] = usedsize;
// 		// Might not always be true, e.g. at edges:
// 		assert(usedsize == ref->data->blocks.beginblock[i]);
// 		usedsize += dims[0] * dims[1] * dims[2];
// 	}

// 	printf("\nthis:\n");
// 	for (int i=0; i<this->data->nrblocks; i++) {
// 		this->get_block_dimensions(i, testdims);
// 		printf("beginblock[%d]: %d\n", i, this->data->blocks.beginblock[i]);
// 		printf("block %d: %d %d %d\n", i, testdims[0], testdims[1], testdims[2]);
// 	}

// 	// // fill beginblock
// 	//  usedsize = 0;
// 	// int dims[3];  // dimensions of the current block
// 	// for (int i=0; i<nrblocks; i++) {
// 	// 	ref->get_block_dimensions(i, dims);
// 	// 	data->blocks.beginblock[i] = usedsize;
// 	// 	usedsize += dims[0] * dims[1] * dims[2];
// 	// }

// 	// set number of blocks
// 	data->nrblocks = nrblocks;

// 	// reallocate tel if necessary
// 	if (usedsize > allocsize) {
// 		allocsize = 2 * usedsize;
// 		//data->blocks->tel = reallocate
// 		reallocate<EL_TYPE>(&(data->blocks.tel), allocsize);
// 		//this->reallocate(allocsize); }
// 	}
// 	if (set_zero) {
// 		// fill tel up with zeros
// 		for (int i=0; i<usedsize; i++) {
// 			data->blocks.tel[i] = 0; }
// 	}
// }


void print_tensorInfo(const struct bookkeeper * keeper,
		const TensorInfo * tensor, int spec)
{
	if (spec == 1 || spec == 0) {
		// print the SiteTensor structure
		fprintf(stdout,"SiteTensor data:\n");
		if (tensor->get_data() != NULL) {
			print_siteTensor(keeper, tensor->get_data());
		}
	}
	
	if (spec == 2 || spec == 0) {
		// print the Symsecs structures
		fprintf(stdout, "Symsectors array:\n");
		struct symsecs * syms[3];
		tensor->get_syms(syms);
		for (int i=0; i<3; i++) {
			if (syms[i] != NULL) {
				print_symsecs(keeper, syms[i], 0);
			}
		}
	}

	if (spec == 3 || spec == 0) {
		fprintf(stdout, "Data info:\n");
		if (tensor->get_data() != NULL) {
			fprintf(stdout, "nr_blocks: %d\n", tensor->get_data()->nrblocks);
		}
		fprintf(stdout, "nr_allocated_blocks: %d\n", tensor->get_nr_allocated_blocks());
		fprintf(stdout, "usedsize: %d\n", tensor->get_usedsize());
		fprintf(stdout, "allocsize: %d\n", tensor->get_allocsize());
	}

	if (spec == 4 || spec == 0) {
		int dims[3];  // dimensions of the current block
		fprintf(stdout, "Block details:\n");
		for (int i=0; i<tensor->get_data()->nrblocks; i++) {
			tensor->get_block_dimensions(i, dims);
			fprintf(stdout, "beginblock[%d]: %d\n", i, tensor->get_data()->blocks.beginblock[i]);
			fprintf(stdout, "block dimensions %d: %d %d %d\n", i, dims[0], dims[1], dims[2]);
			fprintf(stdout, "      symindices %d: %d %d %d\n", i, tensor->get_sym_index(i,0),
					tensor->get_sym_index(i,1), tensor->get_sym_index(i,2));
		}
	}
}


// initialize the memory for a tensorInfo object while
// setting its syms to NULL
void init_dataMemory_tensorInfo(TensorInfo * tensor_info)
{
    // syms will be ented onto external structure
    struct symsecs * symptrs[3] = {NULL, NULL, NULL};
    // siteTensor data is explicitly allocated
    struct siteTensor * data;
    mallocate<struct siteTensor>(&data, 1);
    init_null_siteTensor(data);
    // put everything together in a tensorInfo object
    TensorInfo tensor(data, symptrs, true);
    *tensor_info = tensor;  // == shallow copying
}


void print_TensorInfoPair(const struct bookkeeper * opt_bookie,
		const struct bookkeeper * ref_bookie,
		const struct TensorInfoPair * pair, int spec, char which)
{
	if (which == 'a' || which == 'o') {
		// print the opt TensorInfo
		fprintf(stdout,"TensorInfoPair with\nOPT:\n");
		print_tensorInfo(opt_bookie, &(pair->opt), spec);
	}

	if (which == 'a' || which == 'r') {
		// print the ref TensorInfo
		fprintf(stdout,"REF:\n");
		print_tensorInfo(ref_bookie, &(pair->ref), spec);
	}
}

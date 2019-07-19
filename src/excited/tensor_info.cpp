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


// void TensorInfo::copy_symmetry_layout(const TensorInfo * ref)
// {
// 	// set correct type of tensor
// 	this->is_physical = ref->is_physical;

// 	// set the sectors
// 	for (int i=0; i<3; i++) {
// 		this->syms[i] = ref->syms[i];
// 	}

// 	printf("temporary solution used for symsecs!\n");
// }


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


void TensorInfo::renew_block_layout(bool set_zero)
{
	// number of symmetry sectors in each leg
	int nrSecs[3];
	for (int i=0; i<3; i++) { nrSecs[i] = syms[i]->nrSecs; }
	// number of blocks of the new structure
	int nrblocks = nrSecs[0] * nrSecs[1] * nrSecs[2];
	// set the number of blocks
	data->nrblocks = nrblocks;

	// reallocate qnumbers and beginblock if necessary
	if (nr_allocated_blocks < nrblocks) {
//		reallocate<QN_TYPE>(&(data->qnumbers), nrblocks);
		reallocate<int>(&(data->blocks.beginblock), nrblocks);
		nr_allocated_blocks = nrblocks;
	}

	int counter = 0;
	usedsize = 0;
	// fill beginblocks:
	// -> for all symmetry sectors
	for (int i=0; i < syms[0]->nrSecs; i++) {
		for (int j=0; j < syms[1]->nrSecs; j++) {
			for (int k=0; k < syms[2]->nrSecs; k++)
			{
				// collect the symmetry indices of all legs
				int ids[3] = {i,j,k};
				int dims[3];
				// for all legs
				for (int l=0; l<3; l++) {
					// get the appropriate dimension
					dims[l] = get_block_dimension(l, ids[l]);
				}

				// set the appropriate data
				data->blocks.beginblock[counter] = usedsize;

				// increase the counter and used size
				counter++;
				usedsize += dims[0] * dims[1] * dims[2];
			}
		}
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
// 		//printf("%d -> ", data->qnumbers[i]);
// 		data->qnumbers[i] = translate_qn_address(ref->data->qnumbers[i], ref->syms, this->syms);
// 		assert(data->qnumbers[i] != -1);
// 		//printf("%d due to %d\n", data->qnumbers[i], ref->data->qnumbers[i]);
// 		//data->qnumbers[i] = ref->data->qnumbers[i]; // hier translate qnumber als symsec aangepast!!!
// 	}

// 	// fill beginblock
// 	usedsize = 0;
// 	int dims[3];  // dimensions of the current block
// 	int ids[3];   // symmetry indices of the current block
// 	for (int i=0; i<nrblocks; i++) {
// 		// set the element in beginblock
// 		data->blocks.beginblock[i] = usedsize;

// 		get_sym_indices(i, ids);  //impliciet qnumbers bekijken!!!
// 		// loop over all legs
// 		for (int j=0; j<3; j++) {
// 			// set the appropriate dimension
// 			dims[j] = get_block_dimension(j, ids[j]);
// 		}


// 		// dims[j] = (j == OO_link->leg)? OO_link->OO->get_sdim(ids[j]) :  // warning: dim-convention!
// 		//		get_block_dimension(j, ids[j]); }
// 		// increase the usedsize
// 		usedsize += dims[0] * dims[1] * dims[2];
		

// 		// assert(dims[1] == get_block_dimension(1, ids[1]));
// 		// assert(dims[2] == get_block_dimension(2, ids[2]));
// 		// assert(dims[0] == get_block_dimension(0, ids[0]));
// 		//if (dims[0] != get_block_dimension(0, ids[0])) {printf("WARNING from tensor_info.cpp\n");};


// 		// ERROR: WHEN THE OO CONTAINS (0,0) BLOCKS, THE DIMENSION IS NOT READ OUT CORRECTLY WHEN PRINTING AFTERWARDS!

// 		fprintf(stdout, "beginblock[%d]: %d\n", i, data->blocks.beginblock[i]);
// 		fprintf(stdout, "block dimensions %d: %d %d %d\n", i, dims[0], dims[1], dims[2]);
// 		fprintf(stdout, "      symindices %d: %d %d %d\n", i, ids[0], ids[1], ids[2]);

// 	}

// 	// reallocate tel if necessary
// 	if (usedsize > allocsize) {
// 		allocsize = usedsize;
// 		reallocate<EL_TYPE>(&(data->blocks.tel), allocsize);
// 	}

// 	if (set_zero) {
// 		// fill tel up with zeros
// 		for (int i=0; i<usedsize; i++) {
// 			data->blocks.tel[i] = 0; }
// 	}
// }


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


void print_tensorInfo(const struct bookkeeper ** keepers,
		const TensorInfo * tensor, int spec)
{
	// if (spec == 1 || spec == 0) {
	// 	// print the SiteTensor structure
	// 	fprintf(stdout,"SiteTensor data:\n");
	// 	if (tensor->get_data() != NULL) {
	// 		print_siteTensor(keeper, tensor->get_data());
	// 	}
	// }
	
	if (spec == 2 || spec == 0) {
		// print the Symsecs structures
		fprintf(stdout, "Symsectors array:\n");
		struct symsecs * syms[3];
		tensor->get_syms(syms);
		for (int i=0; i<3; i++) {
			if (syms[i] != NULL) {
				print_symsecs(keepers[i], syms[i], 0);
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

			// THIS CUASES OFF COURSE A SEGMENTATION FAULT BECAUSE THE FUNCTION SHOULD BE ADOPTED TO MORE LOGICAL NUMBERING
			tensor->get_block_dimensions(i, dims);

			fprintf(stdout, "beginblock[%d]: %d\n", i, tensor->get_data()->blocks.beginblock[i]);
			fprintf(stdout, "block dimensions %d: %d %d %d\n", i, dims[0], dims[1], dims[2]);
			fprintf(stdout, "      symindices %d: %d %d %d\n", i, tensor->get_sym_index(i,0),
					tensor->get_sym_index(i,1), tensor->get_sym_index(i,2));
		}
	}
}


void print_tensorInfo(const struct bookkeeper * keeper,
		const TensorInfo * tensor, int spec)
{
	const struct bookkeeper * bookies[3] = {keeper,keeper,keeper};
	print_tensorInfo(bookies, tensor, spec);
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


void print_TensorInfoPair(const struct bookkeeper ** opt_bookies,
		const struct bookkeeper ** ref_bookies,
		const struct TensorInfoPair * pair, int spec, char which)
{
	if (which == 'a' || which == 'o') {
		// print the opt TensorInfo
		fprintf(stdout,"TensorInfoPair with\nOPT:\n");
		print_tensorInfo(opt_bookies, &(pair->opt), spec);
	}

	if (which == 'a' || which == 'r') {
		// print the ref TensorInfo
		fprintf(stdout,"REF:\n");
		print_tensorInfo(ref_bookies, &(pair->ref), spec);
	}
}

void print_TensorInfoPair(const struct bookkeeper * opt_bookie,
		const struct bookkeeper * ref_bookie,
		const struct TensorInfoPair * pair, int spec, char which)
{
	const struct bookkeeper * opt_bookies[3] = {opt_bookie,opt_bookie,opt_bookie};
	const struct bookkeeper * ref_bookies[3] = {ref_bookie,ref_bookie,ref_bookie};
	print_TensorInfoPair(opt_bookies, ref_bookies, pair, spec, which);
}

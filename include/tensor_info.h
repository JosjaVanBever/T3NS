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

#ifndef TENSOR_INFO_H
#define TENSOR_INFO_H


class TensorInfo {
	public:
		// constructor
		// @param: syms is array of pointers to symsecs
		TensorInfo(struct siteTensor * data, struct symsecs ** syms, bool is_psite);

		// getters
		// void get_syms(struct symsecs ** result) const;
		void get_syms(struct symsecs ** result) const;
		// struct symsecs * get_syms() const { return get_sym(0); }
		struct symsecs * get_sym(int leg) const;
		struct siteTensor * get_data() const { return data; }
		int get_usedsize() const { return usedsize; }
		int get_allocsize() const { return allocsize; }
		int get_nr_allocated_blocks() const { return nr_allocated_blocks; }

		// get the bond dimension for a certain leg and symmetry sector
		int get_block_dimension(int leg, int sym_index) const;
		// get the bond dimensions for a certain block and all legs
		void get_block_dimensions(int blocknr, int * dims) const;

		// setters
		void set_sym(struct symsecs * new_sym, int leg) {
			syms[leg] = new_sym; }
		// @ param: 1st element of an array of 3 pointers to symsecs
		void set_syms(struct symsecs ** new_syms) {
			for (int i=0; i<3; i++) { syms[i] = new_syms[i]; }}

		// // renew the symmetry sectors based on the contraction of ref via
		// // the given OO_link
		// // @param:
		// //   ref => contains the uncontracted symsecs
		// //   OO_link => maps the contracted symsec onto the optimizing symsec
		// void renew_symsec_layout(const TensorInfo * ref,
		// 	const struct OverlapObjectLink * OO_link);

		// copy the symmetry sector and tesnsor type from ref
		void copy_symmetry_layout(const TensorInfo * ref);

		// renew all block information; optionally set elements to 0
		// @param:
		//   reference => contains the qnumbers
		//   set_zero = fill the tel array with zeros
		void renew_block_layout(const TensorInfo * reference,
				const struct OverlapObjectLink * OO_link, bool set_zero);

	private:
		// Main data:
		// symmetry structures in order |a>|b><c| (cf. siteTensor)
		struct symsecs * syms[3];  // array of pointers
		// actual siteTensor containing the data
		struct siteTensor * data;
		// is the tensor a physical or branching tensor?
		bool is_physical;

		// Help data:
		// size of tel array used to effectively store elements
		int usedsize;
		// allocated size of the tel array; should be >= usedsize
		int allocsize;
		// allocated number of blocks
		int nr_allocated_blocks;

		// Help functions:
		// get the index of the symmetry sector of a block for a certain leg
		int get_sym_index(int blocknr, int leg) const;
		// get the indices of the symmetry sectors of a block for all legs
		void get_sym_indices(int blocknr, int * result) const;
};


struct TensorInfoPair {
	// tensor of the reference state, e.g. the ground state
	TensorInfo ref;
	// tensor of the state that is being optimized,
	// e.g. the excited state
	TensorInfo opt;
};


void print_tensorInfo(const struct bookkeeper * keeper,
		const TensorInfo * tensor, int specification=0);

void print_TensorInfoPair(const struct bookkeeper * opt_bookie,
		const struct bookkeeper * ref_bookie,
		const struct TensorInfoPair * pair, int specification=0);


#endif
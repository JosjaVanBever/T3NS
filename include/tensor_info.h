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

#include "cpp_interface.h"
#include <stdio.h>


class TensorInfo {
	public:
		// constructor
		TensorInfo(struct siteTensor * data, struct symsecs ** syms);

		// getters
		// void get_syms(struct symsecs ** result) const;
		void get_syms(struct symsecs ** result) const;
		// struct symsecs * get_syms() const { return get_sym(0); }
		struct symsecs * get_sym(int leg) const;
		struct siteTensor * get_data() const { return data; }

		// setters
		void set_sym(struct symsecs * new_sym, int leg) {
			syms[leg] = new_sym; }
		// @ param: 1st element of an array of 3 pointers to symsecs
		void set_syms(struct symsecs ** new_syms) {
			for (int i=0; i<3; i++) { syms[i] = new_syms[i]; }}

		// renew the symmetry sectors based on the contraction of ref via
		// the given OO_link
		// @param:
		//   ref => contains the uncontracted symsecs
		//   OO_link => maps the contracted symsec onto the optimizing symsec
		void renew_symsec_layout(const TensorInfo * ref,
			const struct OverlapObjectLink * OO_link);

		// // renew all block information; optionally set elements to 0
		// // @param:
		// //   reference => contains the qnumbers
		// //   set_zero = fill the tel array with zeros
		// void renew_block_layout(const TensorInfo * reference, bool set_zero)

	private:
		// Main data:
		// symmetry structures in order |a>|b><c| (cf. siteTensor)
		struct symsecs * syms[3];  // array of pointers
		// actual siteTensor containing the data
		struct siteTensor * data;
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
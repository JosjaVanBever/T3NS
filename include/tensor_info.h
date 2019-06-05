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
		TensorInfo(struct siteTensor * data);

		// get the data
		struct siteTensor * get_data() const { return data; }
	private:
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
		const TensorInfo * tensor);

void print_TensorInfoPair(const struct bookkeeper * keeper,
		const struct TensorInfoPair * pair);


#endif
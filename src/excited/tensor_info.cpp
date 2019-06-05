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

TensorInfo::TensorInfo(struct siteTensor * data) : data(data) {}


void print_tensorInfo(const struct bookkeeper * keeper,
		const TensorInfo * tensor)
{
	fprintf(stdout,"SiteTensor data:\n");
	print_siteTensor(keeper, tensor->get_data());
	// fprintf(stdout,"Sumsectors array:\n");
}


void print_TensorInfoPair(const struct bookkeeper * keeper,
		const struct TensorInfoPair * pair)
{
	fprintf(stdout,"TensorInfoPair with\nopt:\n");
	print_tensorInfo(keeper, &(pair->opt));
	fprintf(stdout,"ref:\n");
	print_tensorInfo(keeper, &(pair->ref));
}

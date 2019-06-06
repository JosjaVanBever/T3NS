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

TensorInfo::TensorInfo(struct siteTensor * data, struct symsecs ** syms)
	: data(data)
{
	for (int i=0; i<3; i++) { this->syms[i] = syms[i]; }
}


void TensorInfo::get_syms(struct symsecs ** result) const
{
	for (int i=0; i<3; i++) { result[i] = syms[i]; }
}


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

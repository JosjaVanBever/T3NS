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

// OverlapObject::OverlapObject(struct symsecs * ref, struct symsecs * opt,
// 			int opt_dim) //: ref(ref), opt(opt),
// 			//allocsize(opt_dim)
// {
// 	// // redefine opt_dim if smaller than 0
// 	// if (opt_dim < 0) { 2 * opt->totaldims; }

// 	// // help variable
// 	// int nrBlocks = ref->nrSecs;

// 	// // initialize the sparseblocks
// 	// init_null_sparseblocks(&blocks);
// 	// // -> allocate the tel array
// 	// blocks.tel = (EL_TYPE *) safe_malloc(ref->totaldims * opt_dim, EL_TYPE);
// 	// // -> allocate beginblock
// 	// blocks.beginblock = (int *) safe_malloc(nrBlocks, int);

// 	// // allocate the ldim and sdim array
// 	// ldim = (int *) safe_malloc(nrBlocks, int);
// 	// sdim = (int *) safe_malloc(nrBlocks, int);
// }
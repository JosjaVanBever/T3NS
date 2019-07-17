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

#ifndef CPP_TOOLS
#define CPP_TOOLS

#include "macros.h"

// usage: reallocate<T>(&ptr, size)
template <typename T>
void reallocate(T ** ptr, int size) {
	*ptr = (T *) realloc(*ptr, size * sizeof(T));
}

// usage: mallocate<T>(&ptr, size)
template <typename T>
void mallocate(T ** ptr, int size) {
	*ptr = (T *) safe_malloc(size, T);
}

#endif

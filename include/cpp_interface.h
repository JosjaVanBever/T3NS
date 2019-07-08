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

/*
    @language: both C and C++ are supported
    @file_content: this file contains the interface from an arbitrary
        C program to the OverlapCalculator class that allows for
        outprojection used in the calculation of excited states.
*/


#ifndef CPP_INTERFACE_H
#define CPP_INTERFACE_H


#ifdef __cplusplus
extern "C" {

// declare a C++ interface to get acces to native C functions
// from a C++ context
/*********************START*C++-INTERFACE*********************/

// #include "network.h"
// #include "symsecs.h"
// #include "bookkeeper.h"
    
#include <cstddef>  // declaration of size_t

// native C functions that are called from a C++ context
extern void get_bonds_of_site(int site, int * bonds);
extern bool is_pbond(int bond);
extern void print_symsecinfo(struct symsecs * ss);
extern void print_symsecs(const struct bookkeeper * keeper,
        const struct symsecs *cursymsec, int fci);
extern void bookkeeper_get_symsecs_address(const struct bookkeeper * keeper, 
        struct symsecs ** res, int bond);
extern void bookkeeper_get_symsecs_address_arr(const struct bookkeeper * keeper,
        int n, struct symsecs ** symarr, const int * bonds);
extern void print_siteTensor(const struct bookkeeper * keeper, 
        const struct siteTensor * tens);
extern void init_null_sparseblocks(struct sparseblocks * blocks);
extern void * safe_malloc_helper(long long s, size_t t, const char * typ, 
        const char * file, int line, const char * func);
extern void print_network(void);
extern int get_common_bond(int site1 , int site2);
extern void print_bookkeeper(const struct bookkeeper * keeper, int fci);


/***********************END*C++-INTERFACE*********************/

}
#endif

// declare the C interface to get acces to native C++ functions
// from a C context
/**********************START*C-INTERFACE**********************/

// Structure that contains the internal information of a T3NS
// Several T3NSfill's are used together with their common network
// structure to create OverlapCalculators.
typedef struct T3NSfill {
    // pointer to the first element in an array of siteTensors
    struct siteTensor * data;  // in practice 'const'
    // bookeeper containing the symmetry information of the T3NS
    const struct bookkeeper * bookie;
} T3NSfill;

/************************END*C-INTERFACE**********************/


// declare a OverlapCalculator class or structure
#ifdef __cplusplus
    class OverlapCalculator;
#else
    // declare a C struct to reference the OverlapCalculator
    // class from a C context
    typedef struct OverlapCalculator OverlapCalculator;
#endif


#ifdef __cplusplus
extern "C" {
#endif

// declare a C interface to get acces to native C++ functions
// from a C context
/**********************START*C-INTERFACE**********************/

// Functions declared in overlap_calculator.cpp:

//  -> constructor
extern void init_overlap_calculator(const T3NSfill * opt,
        const T3NSfill * ref, const struct network * netw,
        OverlapCalculator ** result);
//  -> get_result
extern int get_result(OverlapCalculator*);

/************************END*C-INTERFACE**********************/

#ifdef __cplusplus
}
#endif


#endif

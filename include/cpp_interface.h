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

// native C functions that are called from a C++ context
// NO FUNCTIONS YET

/***********************END*C++-INTERFACE*********************/

}
#endif

// declare the C interface to get acces to native C++ functions
// from a C context
/**********************START*C-INTERFACE**********************/

// Structure that contains the internal information of a T3NS
// Several T3NSfill's are used together with their common network
// structure to create OverlapCalculators.
struct T3NSfill {
    struct siteTensor ** data;
    struct bookkeeper * bookie;
};

/************************END*C-INTERFACE**********************/


// declare a OverlapCalculator class or structure
#ifdef __cplusplus
    class OverlapCalculator;
    // contains the implementation of the OverlapCalculator class
    #include "overlap_calculator.h"
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

// C interface to the public functions provided for the
// OverlapCalculator class:
//  -> constructor
// extern void init_overlap_calculator(int test,
//     OverlapCalculator** result);
extern void init_overlap_calculator(struct T3NSfill* opt,
        struct T3NSfill* ref, OverlapCalculator** result);
//  -> get_result
extern int get_result(OverlapCalculator*);

/************************END*C-INTERFACE**********************/

#ifdef __cplusplus
}
#endif


#endif

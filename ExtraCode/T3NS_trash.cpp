// struct siteTensor ** states = malloc((excitation + 1) *
        //         sizeof(struct siteTensor *));
        // states[0] = T3NS[0];
        // states[1] = T3NS[0];
        // states[0].data = T3NS;
        // states[0] = {
        //     .data = &T3NS,
        //     .bookie = &prevbookie
        // };


        // char statefiles[] = 'ground.h5';

        // // Read and continue previous calculation.
        // if (arguments.h5file) {
        //         tic(&chrono, READ_HDF5);
        //         printf(">> Reading %s...\n", arguments.h5file);
        //         if(read_from_disk(arguments.h5file, T3NS, rops)) { return 1; }
        //         minocc = 0;
        //         toc(&chrono, READ_HDF5);
        // }
        
        // while (i < excitation) {
        //     if
        //     if (h5ground)
        //     i++;
        // }
        // toc(&chrono, INIT_OOCALC);

// The overlap calculator does not allocate any of the T3NS data operated on.
// It does allocate its own temporary results.
// class OverlapCalculator;

// #if defined(__STDC__) || defined(__cplusplus)
//   //extern void c_function(Fred*);   /* ANSI C prototypes */
//   extern void init_overlap_calculator(int test, OverlapCalculator** result);
// #else
//   //extern void c_function();        /* K&R style */
//   extern OverlapCalculator** init_overlap_calculator();
// #endif


// declare the C++ interface for native C functions
#ifdef __cplusplus
extern "C" {

/*********************START*C++-INTERFACE*********************/

// native C functions that are called from a C++ context
void implicit_convertion_test();

/***********************END*C++-INTERFACE*********************/

}
#endif

// #ifdef __cplusplus
// extern "C" {
//     void implicit_convertion_test();
// }
// #endif





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


#ifndef CPP_INTERFACE
#define CPP_INTERFACE


/**********************START*C-INTERFACE**********************/

// Structure that contains the internal information of a T3NS
// Several T3NSfill's are used together with their common network
// structure to create OverlapCalculators.
typedef struct T3NSfill {
    struct siteTensor ** data;
    struct bookkeeper * bookie;
} T3NSfill;

/************************END*C-INTERFACE**********************/


// declare the OverlapCalculator class or structure
#ifdef __cplusplus
    class OverlapCalculator;
    // contains the implementation of the OverlapCalculator class
    #include "overlap_calculator.h"
#else
    // declare a C struct to reference the OverlapCalculator
    // class from a C context
    typedef struct OverlapCalculator OverlapCalculator;
#endif


// declare the C interface for native C++ functions
#ifdef __cplusplus
extern "C" {
#endif

/**********************START*C-INTERFACE**********************/

// C interface to the public functions provided for the
// OverlapCalculator class:
//  -> constructor
// extern void init_overlap_calculator(int test,
//     OverlapCalculator** result);
extern void init_overlap_calculator(T3NSfill* opt_t3ns,
    OverlapCalculator** result);
//  -> get_result
extern int get_result(OverlapCalculator*);

/************************END*C-INTERFACE**********************/

#ifdef __cplusplus
}
#endif


// // declare the C++ interface for native C functions
// #ifdef __cplusplus
// extern "C" {

// /*********************START*C++-INTERFACE*********************/

// // native C functions that are called from a C++ context
// void implicit_convertion_test();

// /***********************END*C++-INTERFACE*********************/

// }
// #endif

#ifdef __cplusplus
extern "C" {
    void implicit_convertion_test();
}
#endif


#endif



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


#ifndef CPP_INTERFACE
#define CPP_INTERFACE


#ifdef __cplusplus
extern "C" {

// declare the C++ interface to get acces to native C functions
// from a C++ context
/*********************START*C++-INTERFACE*********************/

// native C functions that are called from a C++ context
void implicit_convertion_test();

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
    struct siteTensor ** data;
    struct bookkeeper * bookie;
} T3NSfill;

/************************END*C-INTERFACE**********************/


// declare the OverlapCalculator class or structure
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

// declare the C interface to get acces to native C++ functions
// from a C context
/**********************START*C-INTERFACE**********************/

// C interface to the public functions provided for the
// OverlapCalculator class:
//  -> constructor
// extern void init_overlap_calculator(int test,
//     OverlapCalculator** result);
extern void init_overlap_calculator(T3NSfill* opt_t3ns,
    OverlapCalculator** result);
//  -> get_result
extern int get_result(OverlapCalculator*);

/************************END*C-INTERFACE**********************/

#ifdef __cplusplus
}
#endif


#endif


//struct symsecs * testsym;

        // bookkeeper_get_symsecs_address(opt_t3ns->bookie, &testsym, CB_nrs[0]);

OverlapCalculator::OverlapCalculator(int test)
{
        this->test = test;
}

:
    opt_bookie(opt_t3ns->bookie), ref_bookie(ref_t3ns->bookie)

    // std::cout << "TensorInfoPair with\nopt:\n";
        // print_siteTensor(keeper, opt.get_data());



    // @TEST
    this->test = opt_t3ns->bookie->nr_bonds;
    std::cout << "Test was set to " << test << std::endl;

    // @TEST
        print_TensorInfoPair(ref_t3ns->bookie, &tensorpairs[i]);

        // @TEST
        std::cout << "tensorpairs[i].opt.data: ";
        print_siteTensor(opt_t3ns->bookie, tensorpairs[i].opt.get_data());

        // @TEST
        print_symsecinfo(COB_syms[0]);
        
        // @TEST
        std::cout << "collected bondnrs: " << CB_nrs[0] << ", " << CB_nrs[1]
                << " and " << CB_nrs[2] << std::endl;

fprintf(stdout, "Nr_bonds was %d\n", states[0].bookie->nr_bonds);
            fprintf(stdout, "Nr_bonds of global is %d\n", bookie.nr_bonds);

            /***************************begin test area************************/
        int CB_nrs[3] = {0, 42, 1};
        struct symsecs COB_syms[3];

        struct symsecs * CTB_syms[3];

        // this function should be modified
        bookkeeper_get_symsecs_arr(&bookie, 3, COB_syms, CB_nrs);
        // bookkeeper_get_symsecs_address(&bookie, &CTB_syms[0], CB_nrs[0]);
        bookkeeper_get_symsecs_address_arr(&bookie, 3, CTB_syms, CB_nrs);

        // CTB_syms[0] = &(bookie.v_symsecs[0]);


        int cheat = bookie.v_symsecs[0].totaldims;
        int cheat2 = COB_syms[0].totaldims;
        fprintf(stdout,"COB before: ");
        print_symsecinfo(&COB_syms[0]);
        fprintf(stdout,"CTB before: ");
        print_symsecinfo(CTB_syms[0]);
        COB_syms[0].totaldims = 20;
        bookie.v_symsecs[0].totaldims = 10;

        fprintf(stdout,"COB after: ");
        // we want to get "bond dimension : 10 (1) in 1 non-empty sectors." here!
        print_symsecinfo(&COB_syms[0]);
        fprintf(stdout,"CTB after: ");
        print_symsecinfo(CTB_syms[0]);

        bookie.v_symsecs[0].totaldims = cheat;
        COB_syms[0].totaldims = cheat2;
        fprintf(stdout,"COB restore: ");
        print_symsecinfo(&COB_syms[0]);
        fprintf(stdout,"CTB restore: ");
        print_symsecinfo(CTB_syms[0]);


        print_siteTensor(&bookie, *(states[0].data));

        /***************************end test area************************/

 fprintf(stdout, "max_dim: %d\n", max_dim);  // @TEST

    symsecMatchers[0].set_test_result(0,4);
    symsecMatchers[1].set_test_result(0,4);
    symsecMatchers[2].set_test_result(4,4);

    const Pair<int> * test = symsecMatchers[0].get_result();
    const Pair<int> * ref1 = symsecMatchers[1].get_result();
    const Pair<int> * ref2 = symsecMatchers[2].get_result();
    for (int i=0; i<4; i++) {
        printf("test: %d %d\n", test[i][0], test[i][1]);
        printf("ref1: %d %d\n", ref1[i][0], ref1[i][1]);
        printf("ref2: %d %d\n", ref2[i][0], ref2[i][1]);
    }
    
    // SymsecMatcher matcher(max_dim);
    // symsecMatchers = {matcher, matcher, matcher};  // matcher is copied

    //            result = (int (*)[2]) safe_malloc(2 * max_size, int); }
    //        const int * get_result()[2] const { return result; };

    // @TEST
        void set_test_result(int startval, int size) {
            this->size = size;
            for (int i=0; i<size; i++)  {
                result[i][0]  = (startval != 0)? startval+i : 0;
                result[i][1]  = (startval != 0)? -startval-i : 0;
            }
        }


        //this->set_syms(&(ref->get_syms()));
        //this->set_sym(contracted_leg, OO_link->OO->get_opt());

        //data->qnumbers = (QN_TYPE*) realloc(data->qnumbers, nrblocks * sizeof(QN_TYPE));
                //data->blocks->beginblock = (int*) realloc(data->blocks->beginblock, nrblocks * sizeof(int));

        // init_null_sparseblocks(&blocks);
        // // -> allocate the tel array
        // blocks.tel = (EL_TYPE *) safe_malloc(ref->totaldims * opt_dim, EL_TYPE);
        
        // // -> allocate beginblock
        // blocks.beginblock = (int *) safe_malloc(nrBlocks + 1, int);


// // help variable
//      int nr_blocks = copy->get_nr_blocks()
//      // sparseblocks
//      deep_copy_sparseblocks(this->blocks, copy->blocks, nr_blocks);
//      // ldim and sdim
//      ldim = (int *) safe_malloc(nr_blocks, int);
//      sdim = (int *) safe_malloc(nr_blocks, int);
//      for (int i=0; i<nr_blocks; i++) {
//              ldim[i] = copy->ldim[i];
//              sdim[i] = copy->sdim[i];
//      }

        // void OverlapObject::swap(OverlapObject& s) noexcept
// {
//      using std::swap;
//      swap(this.mArray,s.mArray);
//      swap(this.mSize ,s.mSize);
// }

           //what if nr blocks change!!!!!!!!!!
                                                                  //cannot be the case!

        // // do extra reallocation just for testing @TEST
        // blocks.beginblock = (int *) realloc(blocks.beginblock, nrBlocks); // @TEST
//      ldim = (int *) realloc(ldim, nrBlocks);// @TEST
//      sdim = (int *) realloc(sdim, nrBlocks);// @TEST

// struct symsecs * symptrs[3];
    // for (int i=0; i<3; i++) {
    //     mallocate<struct symsecs>(&(symptrs[i]), 1);
    //     init_memory_symsecs(&(symptrs[i]));
    // }

         // ldim = (int *) safe_malloc(nr_blocks, int);
        // struct siteTensor * tensor = (struct siteTensor *) safe_malloc(1, struct siteTensor);
        
        // struct siteTensor * data;
        // mallocate<struct siteTensor>(&data, 1);
        // init_null_siteTensor(&data);

        // struct symsecs * symptrs[3];
        // for (int i=0; i<3; i++) {
        //     mallocate<struct symsecs>(&(symptrs[i]), 1);
        //     init_null_symsecs(&(symptrs[i]));
        // }

        // TensorInfo tensor(&data, symptrs, true);
        // tensmem[0] = tensor;

        // struct symsecs * syms;
        // mallocate<struct symsecs>(&syms, 3);
        // for (int i=0; i<3; i++) {
        //     init_null_symsecs(&(syms[i]));
        // }

        
    // }

    // struct symsecs * TEMPsyms[3];
    // struct symsecs helpsyms[3];
    // for (int i=0; i<3; i++) {
    //     TEMPsyms[i] = &(helpsyms[i]);
    //     init_null_symsecs(TEMPsyms[i]);
    // }

    // printf("\nTEST1\n"); fflush(stdout);

    // struct siteTensor TEMPdata;
    // // init_1siteTensor(&TEMPdata, 0, 'n');
    // init_null_siteTensor(&TEMPdata);

    // printf("\nTEST2\n"); fflush(stdout);
    // TensorInfo TEMP(&TEMPdata, TEMPsyms, true);

    // printf("\nTEMP:\n");
    // print_tensorInfo(ref_bookie, &TEMP, 0);




        // get the qnumber and block number 

        //                      QN_TYPE qnum = qntypize_address(ids, syms);  // een counter zou ook volstaan
                                //int block_nr = search_symsec_address(&qnum, syms);

                                

        //                      data->qnumbers[block_nr] = qnum;
        //                      data->beginblock[qntypize_address(ids, syms)] = usedsize;
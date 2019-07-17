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


#include "2site_overlap_calculator.h"


// Declare a C interface to OverlapCalculator subclasses
/**********************START*C-INTERFACE**********************/

int perform_testing_2site(TwoSiteOverlapCalculator* calc)
{
    return calc->perform_testing();
}

// void can_you_find_me(TwoSiteOverlapCalculator* calc)
// {
//     calc->can_you_find_me();
// }

/************************END*C-INTERFACE**********************/


// Set OO equal to the contraction of A and B. 'leg' determines the index
// that is not contracted and which symsec of A will determine OO's blocking.
void TwoSiteOverlapCalculator::set_OO_to_contraction(const TensorInfo * ref,
            const TensorInfo * opt, int open_leg, OverlapObject * OO)
{
    // assert that OO is associated to the open leg
    assert(OO->ref == ref->syms[open_leg] && OO->opt == opt->syms[open_leg]);

    // Indices are permuted (=:per) such that the outermost loop corresponds
    // with the index in the open leg == block index in the OverlapObject
    int per[3] = {open_leg, (open_leg+1)%3,(open_leg+2)%3};
    // Collect the matching irrep pairs for each leg 
    for (int i=0; i<3; i++) {
        symsecMatchers[i].set_matching_symsec_indices(ref, opt, i);
    }

    // modify the block layout of OO
    // HERE IS A SEGMENTATION FAULT GENERATED!
    OO->renew_block_layout(&symsecMatchers[open_leg] , true);

    // // irrep indices for the reference and optimizing blocks currently
    // // contributing
    // int r[3], o[3];

    // // Remark: outer loop is completely parallelizable: no overlap in write
    // //   access since independ blocks are manipulated
    // // For all blocks in OO ...
    // for (int i=0; i<symsecMatchers[open_leg].get_size(); i++) {
    //     r[open_leg] = symsecMatchers[open_leg].get_result()[i][0];
    //     o[open_leg] = symsecMatchers[open_leg].get_result()[i][1];
    //     // and for all irrep indices that contribute to this block
    //     for (int j=0; j<symsecMatchers[per[1]].get_size(); i++) {
    //         r[per[1]] = symsecMatchers[per[1]].get_result()[j][0];
    //         o[per[1]] = symsecMatchers[per[1]].get_result()[j][1];
    //         for (int k=0; k<symsecMatchers[per[2]].get_size(); i++) {
    //             r[per[2]] = symsecMatchers[per[2]].get_result()[k][0];
    //             o[per[2]] = symsecMatchers[per[2]].get_result()[k][1];
    //             // add the given contribution to this block.
    //             double prefactor = get_2leg_contraction_prefactor(
    //                     open_leg, r, ref->get_syms());
    //             add_2leg_contraction(prefactor, open_leg
    //                 ref->get_block_info(r[0],r[1],r[2]),
    //                 opt->get_block_info(o[0],o[1],o[2]),
    //                 OO->get_block_info(r[open_leg]));
    //         }
    //     }
    // }
}


// @ TEST
int TwoSiteOverlapCalculator::perform_testing() {
    // print_network();

    printf("Hurray!\n");

    // printf("\n");
    // for (int i=2; i<3; i++) {
    //     fprintf(stdout, "overlaps[%d]:\n", i);
    //     print_overlap_object(ref_bookie, opt_bookie, &(overlaps[i]));
    // }
    // printf("\n");
/*
    printf("\nBefore:\n");
    for (int i=1; i<2; i++) {
        printf("\n\n%d:", i);
        print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]), 0);
    }

    fprintf(stdout, "\noverlaps[%d]:\n", 1);
    print_overlap_object(ref_bookie, opt_bookie, &(overlaps[2]));
*/
    // tensorpairs[1].opt.set_sym(tensorpairs[0].opt.get_sym(0), 0);
    // overlaps[1].set_opt(tensorpairs[0].opt.get_sym(0));

    // HERE IS A SEGMENTATION FAULT GENERATED!
    set_OO_to_contraction(&(tensorpairs[1].ref), &(tensorpairs[1].opt),
            2, &(overlaps[2]));

    // OverlapObjectLink link;
    // get_internal_link(1, 12, &link);

    // struct symsecs * TEMPsyms[3];
    // struct symsecs helpsyms[3];
    // for (int i=0; i<3; i++) {
    //     TEMPsyms[i] = &(helpsyms[i]);
    //     init_null_symsecs(TEMPsyms[i]);
    // }

    // struct siteTensor TEMPdata;
    // init_null_siteTensor(&TEMPdata);
    // TensorInfo TEMP(&TEMPdata, TEMPsyms, true);

    // printf("\nTEMP:\n");
    // print_tensorInfo(ref_bookie, &TEMP, 0);

    // TEMP.copy_symmetry_layout(&(tensorpairs[1].ref));
    // TEMP.renew_block_layout(&(tensorpairs[1].ref), &link, true);
/*
    printf("\nAfter:\n");
    for (int i=1; i<2; i++) {
        printf("\n\n%d:", i);
        fflush(stdout);
        print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]), 0);
    }

    fprintf(stdout, "overlaps[%d]:\n", 1);
    print_overlap_object(ref_bookie, opt_bookie, &(overlaps[2]));
*/
    // printf("\nTEMP:\n");
    // fflush(stdout);
    // print_tensorInfo(ref_bookie, &TEMP, 0);

    return 0;
}

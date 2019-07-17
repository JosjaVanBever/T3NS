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

#include "overlap_calculator.h"
#include "2site_overlap_calculator.h"

// Declare a C interface to OverlapCalculator subclasses
/**********************START*C-INTERFACE**********************/

void init_overlap_calculator(const T3NSfill * opt, const T3NSfill * ref,
    const struct network * netw, OverlapCalculator ** result, int nsites)
{
    if (nsites == 2) {
        printf("Before new TwoSiteOverlapCalculator\n");
        *result = new TwoSiteOverlapCalculator(opt, ref, netw);
        printf("After new TwoSiteOverlapCalculator\n");
    } else {
        printf("Only 2 site optimizations are for free.\n");
        printf("Extensions to 3 and 4 sites require a doctoral grant.\n");
    }
}

int perform_testing_main(OverlapCalculator* calc)
{
    return calc->perform_testing();
}

// void can_you_find_me_too(OverlapCalculator* calc)
// {
//     calc->can_you_find_me();
// }

/************************END*C-INTERFACE**********************/


// get the maximum bond dimension present in a bookkeeper
int get_max_bond_dimension(const struct bookkeeper * bookie) {
//    print_bookkeeper(bookie, 1);  // @TEST

    // initialize variables
    int result = 0, bond_dim = 0;
    // loop over the virtual bonds
    for (int i = 0; i < bookie->nr_bonds; i++) {
        bond_dim = bookie->v_symsecs[i].totaldims;
//        fprintf(stdout, "bond_dim: %d\n", bond_dim);  // @TEST
        if (result < bond_dim) { result = bond_dim; }
    }
    // loop over the physical bonds
    for (int i = 0; i < bookie->psites; i++) {
        bond_dim = bookie->p_symsecs[i].totaldims;
//        fprintf(stdout, "bond_dim: %d\n", bond_dim);  // @TEST
        if (result < bond_dim) { result = bond_dim; }
    }
//    fprintf(stdout, "result: %d\n", result);  // @TEST
    return result;
}


OverlapCalculator::OverlapCalculator(const T3NSfill * opt_t3ns,
        const T3NSfill * ref_t3ns, const struct network * netw) :
        opt_bookie(opt_t3ns->bookie), ref_bookie(ref_t3ns->bookie),
        network(netw)
{
    // help variables
    // C == Current; O == Optimizing; R == Reference; B == Bond
    // NN = Nearest Neighbouring sites
    int CB_nrs[3];                 // Current Bond numbers
    struct symsecs * COB_syms[3];  // Current Optimizing Bond symsecs
    struct symsecs * CRB_syms[3];  // Current Reference Bond symsecs
    struct symsecs * CO_syms;      // Current optimizing symsecs
    struct symsecs * CR_syms;      // Current reference symsecs

    // do some consistency checks of the input
    assert(opt_t3ns->bookie->nr_bonds == ref_t3ns->bookie->nr_bonds);
    assert(opt_t3ns->bookie->nr_bonds != network->nr_bonds);

    // initialize the tensorpairs
    nr_tensorpairs = network->sites;
    tensorpairs = (struct TensorInfoPair *) safe_malloc(nr_tensorpairs,
            struct TensorInfoPair);
    // fill the tensorpairs
    for (int i=0; i<nr_tensorpairs; i++) {
        // get bond indices
        get_bonds_of_site(i, CB_nrs);

        // @TEST:
        // printf("CB_nrs: %d %d %d\n", CB_nrs[0], CB_nrs[1], CB_nrs[2]);

        // create opt TensorInfo
        bookkeeper_get_symsecs_address_arr(opt_t3ns->bookie, 3, COB_syms, CB_nrs);
        tensorpairs[i].opt = TensorInfo(&(opt_t3ns->data[i]), COB_syms, is_psite(i));
        // create ref TensorInfo
        bookkeeper_get_symsecs_address_arr(ref_t3ns->bookie, 3, CRB_syms, CB_nrs);
        tensorpairs[i].ref = TensorInfo(&(ref_t3ns->data[i]), CRB_syms, is_psite(i));
    }

    // initialize empty overlaps
    nr_overlaps = network->nr_bonds;

    printf("Before overlaps = (OverlapObject *) safe_malloc(nr_overlaps, OverlapObject)\n");

    overlaps = (OverlapObject *) safe_malloc(nr_overlaps, OverlapObject);

    printf("After overlaps = (OverlapObject *) safe_malloc(nr_overlaps, OverlapObject)\n");

    for (int i=0; i<nr_overlaps; i++) {
        bookkeeper_get_symsecs_address(ref_t3ns->bookie, &CR_syms, i);
        bookkeeper_get_symsecs_address(opt_t3ns->bookie, &CO_syms, i);

        printf("Before overlaps[i] = OverlapObject(CR_syms, CO_syms)\n");

        overlaps[i] = OverlapObject(CR_syms, CO_syms);

        printf("After overlaps[i] = OverlapObject(CR_syms, CO_syms)\n");
    }
    overlaps[3].~OverlapObject(); //@TEST

    // look for the maximal bond dimension present in the reference network
    int max_dim = get_max_bond_dimension(ref_t3ns->bookie);

    // initialize the symsecMatchers
    for (int i=0; i<3; i++) {
        symsecMatchers[i].set_size(max_dim);
    }

    // for (int i=0; i<nr_tensorpairs; i++) {
    //     fprintf(stdout,"\n-------------\ntensorpairs[%d]:\n", i);
    //     print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]));
    // }
}


OverlapCalculator::~OverlapCalculator()
{
    safe_free(tensorpairs);
    safe_free(overlaps);
}


// get the OO link attached to who and pointing to other
int OverlapCalculator::get_internal_link(int who, int other,
        OverlapObjectLink * result)
{
    return get_links(who, other, &found_other, NULL, result);
}

// get the OO links attached to who but not pointing to exclude
int OverlapCalculator::get_external_links(int who, int exclude,
        OverlapObjectLink * result)
{
    return get_links(who, exclude, &avoid_other, NULL, result);
}


// is other in (from,to)?
bool found_other(int other, int from, int to)
{
    return from == other || to == other;
}

// is other not in (from, to)?
bool avoid_other(int other, int from, int to)
{
    return from != other && to != other;
}


// Get the appropriate links from who with respect to other,
//   using the criterium check and performing the action question(to, who)
//   to the matches 'to' with the criterium.
int OverlapCalculator::get_links(int who, int other,
        bool (*check)(int,int,int), void (*question)(int,int),
        OverlapObjectLink * result)
{
    // test whether the sites are neighbouring
    if (get_common_bond(who, other) == -1 ) {
        fprintf(stderr, "Warning: sites are not neighbouring.");
    }

    int bond_inds[3];     // bond indices
    int (*bond_link)[2];  // adress of current bond
    int size = 0;         // size of the result

    // get the bond indices for 'who'
    get_bonds_of_site(who, bond_inds);
    // loop over the bond indices
    for (int i=0; i<3; i++) {
        // if the bond is a virtual bond
        // (possibly starting or ending in vacuum)
        if (!is_pbond(bond_inds[i])) {
            // printf("bond_inds[%d]: %d\n", i, bond_inds[i]);
            // ask the tensors linked by the bond
            bond_link = &(network->bonds[bond_inds[i]]);
            int from = (*bond_link)[0];
            // printf("from: %d\n", from);
            int to = (*bond_link)[1];
            // printf("to: %d\n", to);
            // @TEST
            // if (to != -1 && from != -1) {
            //     printf("Eureka!");
            // }
            // if appropriate, add an OO_link to the result
            if (check(other, from, to)) {
                // printf("After check: from: %d\n", from);
                // printf("After check: to: %d\n", to);
                (result[size]).OO = &(overlaps[bond_inds[i]]);
                (result[size]).leg = i;
                // and ask the question if there was one
                if (question) { question(to, who); }
                size++; }}}

    // return the size of the final result
    return size;
}


// @ TEST: common features of OverlapObjectCalculator
int OverlapCalculator::perform_testing() {
    printf("GO BACK TO YOUR HOME COUNTRY!");

    // print_network();

    printf("\n");
    for (int i=2; i<3; i++) {
        fprintf(stdout, "overlaps[%d]:\n", i);
        print_overlap_object(ref_bookie, opt_bookie, &(overlaps[i]));
    }
    printf("\n");

    //printf("data->blocks.beginblock: %d\n", *(tensorpairs[1].opt.get_data()->blocks.beginblock));
    //printf("data->blocks.beginblock: %d\n", tensorpairs[1].opt.get_data()->blocks->beginblock);

    // symsecMatchers[0].set_matching_symsec_indices(&(tensorpairs[4].ref), &(tensorpairs[8].opt), 0);
    // symsecMatchers[1].set_matching_symsec_indices(&(tensorpairs[4].ref), &(tensorpairs[8].opt), 1);
    // symsecMatchers[2].set_matching_symsec_indices(&(tensorpairs[4].ref), &(tensorpairs[8].opt), 2);
    // const Pair<int> * test0 = symsecMatchers[0].get_result();
    // const Pair<int> * test1 = symsecMatchers[1].get_result();
    // const Pair<int> * test2 = symsecMatchers[2].get_result();
    // for (int i=0; i<symsecMatchers[0].get_size(); i++) {
    //     printf("test0: %d %d\n", test0[i][0], test0[i][1]);
    // }
    // for (int i=0; i<symsecMatchers[1].get_size(); i++) {
    //     printf("test1: %d %d\n", test1[i][0], test1[i][1]);
    // }
    // for (int i=0; i<symsecMatchers[2].get_size(); i++) {
    //     printf("test2: %d %d\n", test2[i][0], test2[i][1]);
    // }


    printf("\nBefore:\n");
    for (int i=1; i<2; i++) {
        printf("\n\n%d:", i);
        print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]), 4);
    }

    // tensorpairs[1].opt.set_sym(tensorpairs[0].opt.get_sym(0), 0);
    // overlaps[1].set_opt(tensorpairs[0].opt.get_sym(0));

    OverlapObjectLink link;
    get_internal_link(1, 12, &link);

    struct symsecs * TEMPsyms[3];
    struct symsecs helpsyms[3];
    for (int i=0; i<3; i++) {
        TEMPsyms[i] = &(helpsyms[i]);
        init_null_symsecs(TEMPsyms[i]);
    }

    printf("\nTEST1\n"); fflush(stdout);

    struct siteTensor TEMPdata;
    // init_1siteTensor(&TEMPdata, 0, 'n');
    init_null_siteTensor(&TEMPdata);

    printf("\nTEST2\n"); fflush(stdout);
    TensorInfo TEMP(&TEMPdata, TEMPsyms, true);

    printf("\nTEMP:\n");
    print_tensorInfo(ref_bookie, &TEMP, 0);


    // TEMP.renew_symsec_layout(&(tensorpairs[1].ref), &link);
    TEMP.copy_symmetry_layout(&(tensorpairs[1].ref));

    TEMP.renew_block_layout(&(tensorpairs[1].ref), &link, true);
    fflush(stdout);
    // tensorpairs[2].opt.renew_symsec_layout(&(tensorpairs[1].ref), &link);
    // tensorpairs[2].opt.renew_block_layout(&(tensorpairs[1].ref), true);

    printf("\nAfter:\n");
    for (int i=1; i<2; i++) {
        printf("\n\n%d:", i);
        fflush(stdout);
        print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]), 4);
    }
    printf("\nTEMP:\n");
    fflush(stdout);
    print_tensorInfo(ref_bookie, &TEMP, 0);

    // printf("\nBefore:\n");
    // for (int i=0; i<nr_tensorpairs; i++) {
    //     printf("\n%d:", i);
    //     print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]), 2);
    // }

    // tensorpairs[1].opt.set_sym(tensorpairs[0].opt.get_sym(0), 0);
    // overlaps[1].set_opt(tensorpairs[0].opt.get_sym(0));

    // OverlapObjectLink link;
    // get_internal_link(1, 0, &link);
    // tensorpairs[3].opt.renew_symsec_layout(&(tensorpairs[1].ref), &link);

    // printf("\nAfter:\n");
    // for (int i=0; i<nr_tensorpairs; i++) {
    //     printf("\n%d:", i);
    //     print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]), 2);
    // }

    // print_bookkeeper(ref_bookie, 1);

    // /* initialize random seed: */
    // srand(time(NULL));

    // fprintf(stdout, "Testing link searchers:\n");
    // // int i = rand() % network->sites;
    // // int j = rand() % network->sites;
    // int i = 14;
    // int j = 7;
    // fprintf(stdout, "i is %d and j is %d:\n", i, j);
    // OverlapObjectLink internal[3], external[6];
    // int nr_internal = get_internal_link(i, j, internal);
    // int nr_external = get_external_links(i, j, external);
    

    // fprintf(stdout, "Internal:\n");
    // fprintf(stdout, "nr_internal: %d\n", nr_internal);
    // for (int i=0; i<nr_internal; i++) {
    //     print_overlap_object(ref_bookie, opt_bookie, internal[i].OO);
    // }

    // fprintf(stdout, "External:\n");
    // fprintf(stdout, "nr_external: %d\n", nr_external);
    // for (int i=0; i<nr_external; i++) {
    //     print_overlap_object(ref_bookie, opt_bookie, external[i].OO);
    // }

    // printf("We made it?\n");

    // for (int i=0; i<nr_tensorpairs; i++) {
    //     fprintf(stdout,"\n-------------\ntensorpairs[%d]:\n", i);
    //     print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]));
    // }
    // std::cout << "tensorpairs[i].opt.data: ";
    // print_siteTensor(opt_bookie, tensorpairs[0].opt.get_data());
    return 0;
}

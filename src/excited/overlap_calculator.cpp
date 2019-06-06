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


// Declare a C interface to OverlapCalculator
/**********************START*C-INTERFACE**********************/

void init_overlap_calculator(const T3NSfill * opt, const T3NSfill * ref,
    const struct network * netw, OverlapCalculator ** result)
{
    *result = new OverlapCalculator(opt, ref, netw);
}

int get_result(OverlapCalculator* calc)
{
    return calc->get_result();
}

/************************END*C-INTERFACE**********************/


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
        // create opt TensorInfo
        bookkeeper_get_symsecs_address_arr(opt_t3ns->bookie, 3, COB_syms, CB_nrs);
        tensorpairs[i].opt = TensorInfo(&(opt_t3ns->data[i]), COB_syms);
        // create ref TensorInfo
        bookkeeper_get_symsecs_address_arr(ref_t3ns->bookie, 3, CRB_syms, CB_nrs);
        tensorpairs[i].ref = TensorInfo(&(ref_t3ns->data[i]), CRB_syms);
    }

    for (int i=0; i<nr_tensorpairs; i++) {
        fprintf(stdout,"\n-------------\ntensorpairs[%d]:\n", i);
        print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]));
    }
}


OverlapCalculator::~OverlapCalculator()
{
    safe_free(tensorpairs);
}


// get the OO link attached to who and pointing to other
int OverlapCalculator::get_internal_link(int who, int other,
        OverlapObjectLink * result)
{
    return get_links(who, other, &found_other, NULL, result);
}

// get the OO links attached to who but not pointing to exclude
int OverlapCalculator::get_external_link(int who, int exclude,
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
    int bond_inds[3];     // bond indices
    int (*bond_link)[2];  // adress of current bond
    int size = 0;         // size of the result

    // get the bond indices for 'who'
    get_bonds_of_site(who, bond_inds);
    // loop over the bond indices
    for (int i=0; i<3; i++) {
        // if the bond is a virtual bond
        if (bond_inds[i] != -1) {
            // ask the tensors linked by the bond
            bond_link = &(network->bonds[bond_inds[i]]);
            int from = (*bond_link)[0];
            int to = (*bond_link)[1];
            // if appropriate, add an OO_link to the result
            if (check(other, from, to)) {
                result->OO = &(overlaps[bond_inds[i]]);
                result->leg = i;
                // and ask the question if their was one
                if (question) { question(to, who); }
                size++; }}}

    // return the size of the final result
    return size;
}


// @ TEST
int OverlapCalculator::get_result() {
    for (int i=0; i<nr_tensorpairs; i++) {
        fprintf(stdout,"\n-------------\ntensorpairs[%d]:\n", i);
        print_TensorInfoPair(opt_bookie, ref_bookie, &(tensorpairs[i]));
    }
    // std::cout << "tensorpairs[i].opt.data: ";
    // print_siteTensor(opt_bookie, tensorpairs[0].opt.get_data());
    return 0; };
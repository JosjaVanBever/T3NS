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


// WARNING: links to the symsecs are generated via shallow copies
OverlapCalculator::OverlapCalculator(T3NSfill* opt_t3ns, T3NSfill* ref_t3ns,
            struct network * netw)
{
    // // help types:
    // typedef unordered_map<int, OverlapObjectLink> OOMAP;
    // typedef unordered_map<int, OverlapObjectLink[2]> OO2MAP;
    // help variables
    // C == Current; O == Optimizing; R == Reference; B == Bond
    // NN = Nearest Neighbouring sites
    struct symsecs * COB_syms[3];  // Current Optimizing Bond symsecs
    struct symsecs * CRB_syms[3];  // Current Reference Bond symsecs
    int CB_nrs[3];               // Current Bond numbers
    // struct symsecs * CO_syms;    // Current optimizing symsecs
    // struct symsecs * CR_syms;    // Current reference symsecs
    // struct OverlapObjectLink C_links[3];  // Current OOLinks
    // int CNN_nrs[3];              // Current Nearest Neighbouring sites
    // int CB[2];                   // Current Bond
    // int h1, h2;                  // h = help

    // do some consistency checks of the input
    assert(opt_t3ns->bookie->nr_bonds == ref_t3ns->bookie->nr_bonds);
    assert(opt_t3ns->bookie->nr_bonds != netw->nr_bonds);

    // TEMPORARY
    print_bookkeeper(opt_t3ns->bookie, 0);

    // initialize the tensorpairs
    nr_tensorpairs = netw->sites;
    tensorpairs = (TensorInfoPair *) safe_malloc(nr_tensorpairs, TensorInfoPair);
    // fill the tensorpairs
    for (int i=0; i<nr_tensorpairs; i++) {
        // get bond indices
        get_bonds_of_site(i, CB_nrs);
        // create opt TensorInfo
        bookkeeper_get_symsecs_address_arr(opt_t3ns->bookie, 3, COB_syms, CB_nrs);
        // tensorpairs[i].opt = TensorInfo(opt_t3ns->data[i], COB_syms);
        // create ref TensorInfo
        bookkeeper_get_symsecs_address_arr(ref_t3ns->bookie, 3, CRB_syms, CB_nrs);  //SEG FAULT?!!
        // tensorpairs[i].ref = TensorInfo(ref_t3ns->data[i], CRB_syms);
    }

    // TEMPORARY
    // print_tensorinfopair(&tensorpairs[0], opt_t3ns->bookie);

    // TESTING CODE
    this->test = opt_t3ns->bookie->nr_bonds;
    cout << "Test was set to " << test << endl;
}


OverlapCalculator::~OverlapCalculator()
{
    safe_free(tensorpairs);
}

// OverlapCalculator::OverlapCalculator(int test)
// {
//     this->test = test;
// }
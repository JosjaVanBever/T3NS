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


OverlapCalculator::OverlapCalculator(const T3NSfill * opt_t3ns,
        const T3NSfill * ref_t3ns, const struct network * netw) :
        opt_bookie(opt_t3ns->bookie), ref_bookie(ref_t3ns->bookie)
{
    // help variables
    // C == Current; O == Optimizing; R == Reference; B == Bond
    // NN = Nearest Neighbouring sites
    int CB_nrs[3];                 // Current Bond numbers
    struct symsecs * COB_syms[3];  // Current Optimizing Bond symsecs
    struct symsecs * CRB_syms[3];  // Current Reference Bond symsecs

    // do some consistency checks of the input
    assert(opt_t3ns->bookie->nr_bonds == ref_t3ns->bookie->nr_bonds);
    assert(opt_t3ns->bookie->nr_bonds != netw->nr_bonds);

    // initialize the tensorpairs
    nr_tensorpairs = netw->sites;
    tensorpairs = (struct TensorInfoPair *) malloc(nr_tensorpairs *
            sizeof(struct TensorInfoPair));
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
    free (tensorpairs);
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
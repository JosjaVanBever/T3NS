// TRASH

// // i is the index in the corresponding ref->blocks
// OverlapObject::get_index(int i, int leg) {
// 	return indexize()[leg]
// 		{a0,a1,a2},A->tensor->qnnumbers[a0],
// 						&(A->tensor->symsec[leg]))[leg];
// }

// //////////////////////////////////////////////////////////////////////

// index = indexize({a0,a1,a2},A->tensor->qnnumbers[a0],
// 						&(A->tensor->symsec[leg]));


// // tensor A is considered as the reference for OO's symsectioning
// void set_OverlapObject(const TensorInfo * A, const TensorInfo * B, int leg,
// 			OverlapObject * OO) {
// 	OO->reset_blocks();
// 	// the looping indices are permuted such that the first one corresponds with 'leg'
// 	int per[3] = {leg%3,(leg + 1)%3,(leg + 2)%3}; // per == permutation
// 	for (int i=0; i<3; i++) {
// 		symsecMatchers[per[i]].set_matching_symsec_indices();
// 	}
// 	for (int i=0; i<symsecMatchers[per[0]].get_size(); i++) {
// 		a0 = symsecMatchers[per[0]].get_result()[i][0];
// 		b0 = symsecMatchers[per[0]].get_result()[i][1];
// 		for (int j=0; j<symsecMatchers[per[1]].get_size(); i++) {
// 			a1 = symsecMatchers[per[1]].get_result()[j][0];
// 			b1 = symsecMatchers[per[1]].get_result()[j][1];
// 			for (int k=0; k<symsecMatchers[per[2]].get_size(); i++) {
// 				a2 = symsecMatchers[per[2]].get_result()[k][0];
// 				b2 = symsecMatchers[per[2]].get_result()[k][1];
// 				// do the actual block manipulation
// 				index = indexize({a0,a1,a2},A->tensor->qnnumbers[a0],
// 						&(A->tensor->symsec[leg]));
// 				OO->add_product_to_block(get_prefactor(i), A->get_block(a0,a1,a2),
// 						B->get_block(b0,b1,b2), a0);
// 			}
// 		}
// 	}
// }



// 				// do the actual block manipulation
// 				OO->add_product_to_block(get_prefactor(),
// 						A->get_block_info(a[0],a[1],a[2]),
// 						B->get_block_info(b[0],b[1],b[2]), a[leg]);

// inline void OverlapObject::add_product_to_block(double prefactor,
// 		const struct BlockInfo A, const struct BlockInfo B; int index) {
// 	struct BlockInfo OO = this->get_block_info(index);
// 	struct contractinfo cinfo = {
// 		.tensneeded = {0, 1, 2},
// 		. = {CblasNoTrans, CblasNoTrans}
// 		// TODO
// 		.M = A.,
// 		.N = B.,
// 		.K = A.,
// 		.L = ,
// 		.lda = ,
// 		.ldb = ,
// 		.ldc = ,
// 		.stride = ,
// 	};
// 	EL_TYPE tels[3] = {OO->start, A->start, B->start};
// 	do_contract(cinfo, tels, prefactor, 1);
// }


// result->blocks.tel["begin"block[i]] = prefactor * data->blocks.tel["begin"block[i] *
// 				toAdd->OO->blocks.tel["begin"block[j]];

// // TODO
// 		result->blocks.tel["begin"block[i]] = prefactor * data->blocks.tel["begin"block[i] *
// 				toAdd->OO->blocks.tel["begin"block[j]];


// initialize the internal_links and external_links
	// internal_links.reserve(nr_tensorpairs, nr_tensorpairs);
	// external_links.reserve(nr_tensorpairs, nr_tensorpairs);
	// 	// REMARK: this is really inefficient, but for x := nr_tensorspairs
	// 	// = O(100), a memory requirement of O(x²) is not such a problem
	// for (int i=0; i< nr_tensorpairs; i++) {
	// 	for (int j=0; j< nr_tensorpairs; j++) {
	// 		internal_links(i,j) = null;
	// 		external_links(i,j) = null;
	// 	}
	// 	int bond_nrs[3];
	// 	get_bonds_of_site(i, bond_nrs);
	// 	for (int k=0; k<2; k++) {
	// 		int * CB = ref_t3ns->bonds[k]  // CB = current_bond
	// 		int other = (CB[0] != i)? CB[0] : CB[1];
	// 		struct OverlapObjectLink current_link = {
	// 			.OO = &overlaps[other],
	// 			.leg = k
	// 		};
	// 		if (other == j) {
	// 			internal_links[i][j] = current_link;
	// 		} else {
	// 			external_links[i][j] = current_link;
	// 			}
	// 		}
	// 	}


// for (int a=0; a<2; a++) {
// 			// get the current bond number
// 			int CB_nr = CB_nrs[a];
// 			// skip when the current bond is a physical bond
// 			if (! is_pbond(CB_nr)) {
// 				// get the bond for the current bond index
// 				int * CB = ref_t3ns->bonds[CB_nr];
// 				// look up the number of the other tensor at this bond
// 				int other_nr = (CB[0] != i)? CB[0] : CB[1];
// 				// create a link to the OO associated to the current bond index
// 				struct OverlapObjectLink current_link = {
// 					.OO = &overlaps[CB_nr],
// 					.leg = a  // order convention of get_bonds_of_site
// 				}; // WRITING
// 				for (int j=0; j<2; j++)
// 				if (other_nr == nn_nrs[j]) {
// 					internal_links[i][j] = current_link;
// 				} else {
// 					external_links[i][j] = current_link;
// 					}
// 				}
// 			}

// }
// 	internal_links

	// internal_links = safe_malloc(nr_tensorpairs, OverlapObjectLink *);
	// for (int i=0; i< nr_tensorpairs; i++) {
	// 	internal_links = safe_malloc(nr_tensorpairs, OverlapObjectLink);
	// }
	
	// for (int b=0; b<3; b++) {
			
	// 		// current_bonds[b] = &(opt_t3ns->bonds[current_bond_nrs[b]])
	// 	}

/*
template<class T>
class NarrowDataBlock : public DataBlock {
	public:
		// d1 == large dimension; d2 == small dimension
		NarrowDataBlock(size_t dlarge=0, size_t dsmall=0) : d1(dlarge), d2(dsmall) {
		data = safe_malloc(dlarge*dsmall, pair<int,T>); }
		// access operators
        // @param: i,j and k are the indexes where the element is stored
        T * operator()(int i, int j) { 
        	for (int h=i*d1; h < i*d1+d2; h++ ) {
        		if (data[h].key() == j) {
        			return data[h].value();
        		}
        	}

}
*/

template <class T>
class DataBlock {
	public:
		// con-/destructor
		DataBlock(size_t d1=0, size_t d2=0) : d1(d1), d2(d2) {
			data = safe_malloc(d1*d2, T); }
		void reserve(size_t d1, size_t d2) {
			this->d1 = d1;
			this->d2 = d2;
			data = safe_malloc(d1*d2, T); }
		~DataBlock() { free(data); }

		// getters for the dimensions adn data
		size_t get_d1() const { return d1; }
        size_t get_d2() const { return d2; }
        T * get_data() { return data; }

        // access operators
        // @param: i,j and k are the indexes where the element is stored
        T * operator()(int i, int j) { return data[i + j*d2]; }
        T const * operator()(int i, int j) const { return data[i + j*d2]; }

    private:
    	size_t d1,d2;  // the dimensions of the data block
    	T * data;      // the data stored in the data block
}

class sth {
	std::unordered_map<int,int> data;
	data[key] = value
	// == mini dictionary

}
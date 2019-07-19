// calculate all OO's so that they point towards optimization_center
void TwoSiteOverlapCalculator::prepare_first_run(int * optimization_center)
{
	// for each of the tensor pairs involved
	for (int i=0; i<2; i++) {
		// prepare its OOs
		prepare_OOs_for_pair(optimizing_tensors[i], optimizing_tensors[1-i]);
		// remember its index
		last_optimized[i] = optimizing_tensors[i];
	}
}


// calculate all OO's at the bonds of the tensor pair who excuding the one
// at the bond between who and exclude
prepare_OOs_for_pair (int who, int exclude)
{
	// If we're living on the edge:
	if (who == -1) {

		// get the OO living here
		OverlapObjectLink to_calc;
		get_internal_link(who, exclude, to_calc);

		// get the pair at 'exclude'
		TensorInfoPair * pair = tensorpairs[exclude];

		// The OO is simply the contraction of pair->ref and pair->opt.
		set_OO_to_contraction(pair.get_ref(), pair.get_opt(), to_calc->leg,
				to_calc->OO);
	
	// Else if there are some neighbours to contact first:
	} else {

		// collect OO links to them while asking to update their OO's
		OverlapObjectLink[2] external;
		size = get_links(who, exclude, &avoid_other,
				&prepare_OOs_for_pair, external);

		// get the pair at 'who'
		TensorInfoPair * pair = tensorpairs[who];

		// use the links to merge their OO's with the reference tensor
		contract_reference_with_OO_list(pair.get_ref(), external,
				&(tensbmem[0]), &(tensmem[0]));

		// merge this result with the optimizing tensor at who
		set_OO_to_contraction(&(tensmem[0]), pair.get_opt(), to_calc->leg,
				to_calc->OO);
	}
}

// WARNING: WHAT IF SWITCHING DIRECTION?

// Compare the current (=:C) and previous (=:P) optimization center and update
// the OO that is directed from the tensor (in P but not in C) to the tensor
// (in P and in C).
void TwoSiteOverlapCalculator::prepare_OOs(int * optimizing_tensors) {
	int * common = null;      // common index with P
	int * left_over = null;   // index that's only in P
	TensorInfo * memory;      // memory location of intermediate saved result

	// compare the current and previous optimization center
	for (int i=0; i<2; i++) {
		int current = last_optimized[i];
		if(linSearch(current, optimizing_tensors, 2, sort_int[1],
					sizeof(int)) == -1) {
			memory = &(tensmem[i]);
			*left_over = current;
		} else {
			*common = current; }}

	// A necessary condition is that the current and previous optimization
	// center have at least one tensor in common.
	assert(common != null);

	// Update if the optimization center has moved
	if (left_over != null) {
		update(*left_over, *common, memory); }
}


// The optimizing tensors are assumed to correspond to the indices of the
// network used in the constructor. No further assumptions are made.
TwoSiteOverlapCalculator::get_overlap_vector(int * optimizing_tensors,
			struct siteTensor * result)
{
	// Firstly, the Overlap Object between the previously updated sites might
	// be updated so that it can be used in the calculation of the new overlap
	// vector. This means that the particular OO needs to be directed to the
	// current optimization center.
	prepare_OOs(optimizing_tensors);

	// Secondly, for the tensorpair at tensorpairs[optimization_tensors[i]] the
	// contraction op its reference tensor with the OO's at the border of the
	// optimization center is stored in tensmem[i].

	// for each of the tensors involved
	for (int i=0; i<2; i++) {
		// classify the 2 tensor indices involved
		int current_index = optimizing_tensors[i]
		int other_index = optimizing_tensors[1-i]
		// get the appropriate reference site and overlap object
		TensorInfo * reference = tensorpairs[current_index].get_ref();
		OverlapObjectLink external[2] = external_links[current_index,other_index];
		// perform the actual contraction
		contract_reference_with_OO_list(reference, external, &(tensbmem[i]), &(tensmem[i]));
	}

	// Thirdly, the contractions in tensmem[0] and tensmem[1] are put together
	// in a two site tensor. Considering the indices of this object as one index,
	// this is the overlap vector needed.
	auto T3NS = { tensmem[0].get_data(), tensmem[1].get_data(), -1};
	makesiteTensor(result , T3NS, optimizing_tensors, 2)

	// Finally, remember the optimization center used for subsequent calculations.
	for (int i=0; i<2; i++) { last_optimized[i] = optimizing_tensors[i]; }
}


void OverlapCalculator::contract_reference_with_OO(const TensorInfo * reference,
			const OverlapObjectLink * OO_link, TensorInfo * result)
{
	...

	// perform the actual block contractions
	int contracted_leg = OO_link->leg;
	for (int i=0; i < reference->data->nrblocks; i++) {
		double prefactor = get_1leg_contraction_prefactor(); 
		// labeling in the OO is done using a single leg 
		j = tensor->get_leg_index(i,contracted_leg);
		add_1leg_contraction(prefactor, contracted_leg,
				reference->get_block_info(i),
				OO_link->OO->get_block_info(j),
				result->get_block_info(i));
	}
}


void OverlapCalculator::contract_reference_with_OO_list(const TensorInfo * ref,
			const OverlapObjectLink * OO_link_list, TensorInfo * memory,
			TensorInfo * result)
{
	if (ref->is_physical_tensor()) {
		contract_reference_with_OO(ref, OO_link_list[0], result);
	} else {
		contract_reference_with_OO(ref, OO_link_list[0], memory);
		contract_reference_with_OO(memory, OO_link_list[1], result);
	}
}


// Update the OO present at the bond between the tensorpairs at index from
// and to. The OO is updated such that it is directed towards 'to'.
// Memory contains the contraction of from->ref and its neighbouring OO.
void TwoSiteOverlapCalculator::update(int from, int to, const TensorInfo * memory) {
	// get the tensorinfo's and OO link for the current update
	TensorInfo * opt = &(tensorpairs[from].opt);
	OverlapObjectLink * internal = internal_link[from, to];
	// Save in internal->OO the contraction of memory and the relevant
	// site of the optimizing T3NS.
	set_OO_to_contraction(memory, opt, internal->leg, internal->OO);
}


// Set OO equal to the contraction of A and B. 'leg' determines the index
// that is not contracted and which symsec of A will determine OO's blocking.
void TwoSiteOverlapCalculator::set_OO_to_contraction(const TensorInfo * ref,
			const TensorInfo * opt, int open_leg, OverlapObject * OO)
{
	...

	// Remark: outer loop is completely parallelizable: no overlap in write
	//   access since independ blocks are manipulated
	// For all blocks in OO ...
	for (int i=0; i<symsecMatchers[open_leg].get_size(); i++) {
		r[open_leg] = symsecMatchers[open_leg].get_result()[i][0];
		o[open_leg] = symsecMatchers[open_leg].get_result()[i][1];
		// and for all irrep indices that contribute to this block
		for (int j=0; j<symsecMatchers[per[1]].get_size(); i++) {
			r[per[1]] = symsecMatchers[per[1]].get_result()[j][0];
			o[per[1]] = symsecMatchers[per[1]].get_result()[j][1];
			for (int k=0; k<symsecMatchers[per[2]].get_size(); i++) {
				r[per[2]] = symsecMatchers[per[2]].get_result()[k][0];
				o[per[2]] = symsecMatchers[per[2]].get_result()[k][1];
				// add the given contribution to this block.
				double prefactor = get_2leg_contraction_prefactor(
						open_leg, r, ref->get_syms());
				add_2leg_contraction(prefactor, open_leg
					ref->get_block_info(r[0],r[1],r[2]),
					opt->get_block_info(o[0],o[1],o[2]),
					OO->get_block_info(r[open_leg]));
			}
		}
	}
}

inline struct BlockInfo TensorInfo::get_block_info(int index) {
	int dims[3];
	get_block_dimensions(index, dims);
	return struct BlockInfo {
		.start = data->blocks.beginblock[index],
		.d1 = dims[0],
		.d2 = dims[2],
		.d3 = dims[3]
	}
}


inline struct BlockInfo TensorInfo::get_block_info(int i, int j, int k) {
	int index = search_symsec(qntypize({i,j,k}, syms));
	return struct BlockInfo {
		.start = data->blocks.beginblock[index],
		.d1 = get_block_dimension(i),
		.d2 = get_block_dimension(j),
		.d3 = get_block_dimension(k)
	}
	return get_block_info(index);
}



/*********************************Varia**********************************/

// information about a block that starts at memory location 'start' and
// has dimensions d1, d2 and d3. If the block is column major, the index
// of element (i,j,k) is given by i + d1 * j + d1 * d2 * k
typedef struct BlockInfo {
	EL_TYPE * start;  // pointer to first element
	int d1, d2, d3;   // dimensions of the block
} BlockInfo;


double get_1leg_contraction_prefactor()
{
	return 1;
}


// All prefactors are trivial if one uses a gauge in which all tensors
// upstream the renormalization flow are right canonical and all others
// are left canonical. (move around hollow and filled dots to achieve this)
double get_2leg_contraction_prefactor()
{
	return 1;
}


// Map a value onto -1 if odd and onto +1 if even.
// This is equivalent to (-1)**(value).
int grade(int value)
{
	return (value % 2 == 0)? 1 : -1;
}
int grade(double value)
{
	return (int(round(value)) % 2 == 0)? 1 : -1;
}


// returns the value of an irrep for a given symmetry and symmetry sector
int get_irrep_value(enum symmetrygroup sym, int irrep_index,
		const struct symsecs * sectors, const struct bookkeeper * bookie)
{
	for (int i = 0; i < bookie.nrSyms; i++) {
	    if (bookie.sgs[i] == sym) {
	        return sectors->irreps[irrep_index][i]; }}

	// Symmetry is not present in the bookkeeper.
	fprintf(stderr, "Warning: no %s is used.\n", sym);
	return -1;
}


// returns 0 (even parity) or 1 (odd parity)
int get_parity(int irrep_index, const struct symsecs * sectors,
		const struct bookkeeper * bookie)
{
	return get_irrep_value(Z2, irrep_index, sector, bookie);
}


// returns 2j
int get_angular_momentum(int irrep_index, const struct symsecs * sectors,
		const struct bookkeeper * bookie)
{
	return get_irrep_value(SU2, irrep_index, sector, bookie);
}

// WARNING: TO MUCH FREEDOM FOR DIMENSIONS
// Perform one of the following operations, taking into account 
// possible mismatch in dimensions:
//   if open_leg == 0:
//     C_{ij} += prefactor * A_{i(kl)} * B_{j(kl)}
//   if open_leg == 1:
//     C_{ij} += prefactor * A_{kil} * B_{kjl}
//   if open_leg == 2:
//     C_{ij} += prefactor * A_{(kl)i} * B_{(kl)j}
void add_2leg_contraction(double prefactor, int open_leg,
		const BlockInfo * A, const BlockInfo * B, BlockInfo * C)
{
	// To take into account mismatch in dimensions, it suffices
	// to take for K and L the minimal contraction dimension available.
	if (open_leg == 0) {
		// C_{ij} += prefactor * A_{i(kl)} * B_{j(kl)}
		int K = min(A->d2 * A->d3, B->d2 * B->d3);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
				A->d1, B->d1, K, prefactor, A->start, A->d1,
				B->start, B->d1, 1, C->start, C->d1);
	}
	else if (open_leg == 1) {
		// C_{ij} += prefactor * A_{kil} * B_{kjl}
		EL_TYPE (*tensors)[3] = {A->start, B->start, C->start};
		struct contractinfo cinfo = {
			.tensneeded = {0, 1, 2},
			.trans = {CblasTrans, CblasNoTrans},
			.M = A->d2,
			.N = B->d1,
			.K = min(A->d1, B->d1),
			.L = min(A->d3, B->d3),
			.lda = A->d1,
			.ldb = B->d1,
			.ldc = C->d1,
			.stride = {A->d1 * A->d2, B->d1 * B->d2 , 0}
		}
		do_contract(&cinfo, tensors, prefactor, 1);
	}
	else if (open_leg == 2) {
		// C_{ij} += prefactor * A_{(kl)i} * B_{(kl)j}
		int K = min(A->d1 * A->d2, B->d1 * B->d2);
		cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
				A->d3, B->d3, K, prefactor, A->start, A->d1 * A->d2,
				B->start, B->d1 * B->d2, 1, C->start, C->d1);
	}
}


// WARNING: TO MUCH FREEDOM FOR DIMENSIONS
// Perform one of the following operations, taking into account 
// possible mismatch in dimensions:
//   if contracted_leg == 0:
//     C_{ijk} += prefactor * B_{li} * A_{ljk}
//   if contracted_leg == 1:
//     C_{ijk} += prefactor * A_{ilk} * B_{lj}
//   if contracted_leg == 2:
//     C_{ijk} += prefactor * A_{ijl} * B_{lk}
void add_1leg_contraction(double prefactor, int contracted_leg,
		const BlockInfo * A, const BlockInfo * B, BlockInfo * C)
{
	if (contracted_leg == 0) {
		// C_{ijk} += prefactor * B_{li} * A_{ljk}
		int K = min(A->d1, B->d1);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
				B->d2, A->d2 * A->d3, K, prefactor, A->start, A->d1,
				B->start, B->d1, 1, C->start, C->d1);
	}
	else if (contracted_leg == 1) {
		// C_{ijk} += prefactor * A_{ilk} * B_{lj}
		EL_TYPE (*tensors)[3] = {A->start, B->start, C->start};
		struct contractinfo cinfo = {
			.tensneeded = {0, 1, 2},
			.trans = {CblasNoTrans, CblasNoTrans},
			.M = A->d1,
			.N = B->d2,
			.K = min(A->d2, B->d1),
			.L = A->d3,
			.lda = A->d1,
			.ldb = B->d1,
			.ldc = C->d1,
			// WARNING: CHECK CORRECTION AFTER AFIRMATION
			.stride = {A->d1 * A->d2, 0, B->d1 * B->d2}
		}
		// dimensional checks
		assert(B->d2 == C->d2);
		assert(A->d3 == C->d3);
		do_contract(&cinfo, tensors, prefactor, 1);
	}
	else if (contracted_leg == 2) {
		// C_{ijk} += prefactor * A_{ijl} * B_{lk}
		int K = min(A->d3, B->d1);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
				A->d1 * A->d2, B->d2, K, prefactor, A->start,
				A->d1 * A->d2, B->start, B->d1 * B->d2, 1,
				C->start, C->d1 * C->d2);
	}
}


// get the info for a block contraction
inline struct BlockInfo OverlapObject::get_block_info(int index)
{
	return struct BlockInfo {
		.start = this->blocks.beginblock[index],
		.d1 = ldim[index],
		.d2 = sdim[index],
		.d3 = 0
	}
}
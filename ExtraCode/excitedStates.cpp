
#include <stdexcept>

/*********************************OverlapCalculator**************************/

typedef struct T3NSfill {
	siteTensor ** data;
	struct bookkeeper * bookie;
} T3NSInfo;

#include <unordered_map>


// The overlap calculator does not allocate any of the T3NS data operated on.
// It does allocate its own temporary results.
class OverlapCalculator {
	public:
		// constructor and destructor
		OverlapCalculator(T3NSfill* opt, T3NSfill* ref, struct network * netw);
		~OverlapCalculator();
		// get the overlap vector used for the Heff diagonalization
		virtual void get_overlap_vector(int * optimization_center, double * result);
		// update the internal data structure before the next diagonalization
		virtual void update();
	private:
		// main data
		TensorInfoPair * tensorpairs;
		int nr_tensorpairs;  // length of tensorpairs
		OverlapObject * overlaps;
		int nr_overlaps;  // length of overlaps
		// help structures to search relevant OverlapObjects:
		// ->  internal_links[x][y] == link to the OO associated to the bond between
		//     tensorpairs[x] and tensorpairs[y], and with reference to tensorpairs[x]
		unordered_map<int, OverlapObjectLink> * internal_links;
		// ->  external_links[x][y] == links to the OO's "before" tensorpairs[x] (=:tx)
		//     with respect to tensorpairs[y] (=: ty). I.e. the OO's used to update
		//     internal_links[x][y] when the sweep direction goes from tx to ty.
		unordered_map<int, OverlapObjectLink[2]> * external_links;
		// help structures to search common irreps
		SymsecMatcher[3] symsecMatchers;
};

// WARNING: links to the symsecs are generated via shallow copies
OverlapCalculator::OverlapCalculator(T3NSfill* opt_t3ns, T3NSfill* ref_t3ns,
			struct network * netw) :
{
	// help types:
	typedef unordered_map<int, OverlapObjectLink> OOMAP;
	typedef unordered_map<int, OverlapObjectLink[2]> OO2MAP;
	// help variables
	// C == Current; O == Optimizing; R == Reference; B == Bond
	// NN = Nearest Neighbouring sites
	struct symsecs COB_syms[3];  // Current Optimizing Bond symsecs
	struct symsecs CRB_syms[3];  // Current Reference Bond symsecs
	int CB_nrs[3];               // Current Bond numbers
	struct symsecs * CO_syms;    // Current optimizing symsecs
	struct symsecs * CR_syms;    // Current reference symsecs
	struct OverlapObjectLink C_links[3];  // Current OOLinks
	int CNN_nrs[3];              // Current Nearest Neighbouring sites
	int CB[2];                   // Current Bond
	int h1, h2;                  // h = help

	// do some consistency checks of the input
	assert(opt_t3ns->bookie->nr_bonds == ref_t3ns->bookie->nr_bonds);
	assert(opt_t3ns->bookie->nr_bonds != netw->nr_bonds);

	// initialize the tensorpairs
	nr_tensorpairs = netw->sites;
	tensorpairs = safe_malloc(nr_tensorpairs, TensorInfoPair);
	// fill the tensorpairs
	for (int i=0; i<nr_tensorpairs; i++) {
		// get bond indices
		get_bonds_of_site(i, CB_nrs);
		// create opt TensorInfo
		bookkeeper_get_symsecs_arr(opt_t3ns->bookie, 3, COB_syms, CB_nrs);
		tensorpairs[i].opt = TensorInfo(opt_t3ns->data[i], COB_syms);
		// create ref TensorInfo
		bookkeeper_get_symsecs_arr(ref_t3ns->bookie, 3, CRB_syms, CB_nrs);
		tensorpairs[i].ref = TensorInfo(ref_t3ns->data[i], CRB_syms);
	}

	// initialize empty overlaps
	nr_overlaps = netw->nr_bonds;
	overlaps = safe_malloc(nr_overlaps, OverlapObject);
	for (int i=0; i<nr_overlaps; i++) {
		bookkeeper_get_symsecs(ref_t3ns->bookie, CR_syms, i);
		bookkeeper_get_symsecs(ref_t3ns->bookie, CO_syms, i);
		overlaps[i] = OverlapObject(CR_syms, CO_syms);
	}

	// initialize the symsecMatchers
	// TODO: extract maximal size max_size needed from reference T3NS
	SymsecMatcher matcher(max_size);
	symsecMatchers = {matcher, matcher, matcher};  // matcher is copied

	// initialize the internal_links and external_links
	internal_links = safe_malloc(nr_tensorpairs, OOMAP);
	external_links = safe_malloc(nr_tensorpairs, OO2MAP);
	// fill the internal_links and external_links
	for (int i=0; i< nr_tensorpairs; i++) {
		// get bond indices
		get_bonds_of_site(i, CB_nrs);
		// get the indices of the corresponding tensor neighbours
		// and the OOLinks to the corresponding OverlapObjects
		for (int a=0; a<2; a++) { 
			// if the current bond number points to a virtual bond
			if (! is_pbond(CB_nr)) {
				// get the bond for the current bond index
				CB = netw->bonds[CB_nrs[a]];
				// look up the index of the other tensor at this bond
				CNN_nrs[a] = (CB[0] != i)? CB[0] : CB[1];
				// store an OOLink to the OO associated to the current bond
				C_links[a] = {
				.OO = &overlaps[CB_nrs[a]],
				.leg = a  // order convention of get_bonds_of_site
			} else {
				// set to default values
				CNN_nrs[a] = -1;
				C_links[a] = null;
			}
		}
		// fill the internal_links and external_links for the current tensor
		for (int a=0; a<2; a++) {
			// fill the internal link
			internal_links[i][CNN_nrs[a]] = C_links[a];
			// take care that only the last element of the external link can be null
			h1 = (C_links[(a+1)%3] != null)? 1 : 2; h2 = 1 - h1;
			// fill the external link
			external_links[i][CNN_nrs[a]] = {C_links[(a+h1)%3], C_links[(a+h2)%3]};
		}
	}
}


OverlapCalculator::~OverlapCalculator() {
	free(tensorpairs); free(overlaps); free(internal_links); free(external_links);
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
	return get_links(who, other, &avoid_other, NULL, result);
}


// is other in (from,to)?
bool found_other(other, from, to)
{
	return from == other || to == other;
}


// is other not in (from, to)?
bool avoid_other(other, from, to)
{
	return from != other && to != other;
}


// Get the appropriate links from who with respect to other,
// using the criterium check and performing the action question(to, who)
// to the matches 'to' with the criterium.
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
		if (NONSENS bond_inds[i] != -1) {
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



class TwoSiteOverlapCalculator : public OverlapCalculator {
	public:
		// TODO constructor and destructor
		TwoSiteOverlapCalculator(T3NSInfo* opt, T3NSInfo* ref,
				struct network * netw) : OverlapCalculator(opt,ref,netw),
				last_optimized({-1,-1}) {
					for (int i=0; i<2; i++) {
						init_null_siteTensor(tensmem[i]);
						init_null_siteTensor(tensbmem[i]);
					}
				};
		// get the overlap vector used for the Heff diagonalization
		void get_overlap_vector(int * optimization_center, double * result);
		// update the internal data structure before the next diagonalization
		void update();
	private:
		void set_OO_to_contraction(const TensorInfo * A, const TensorInfo * B,
				int leg, OverlapObject * OO);
		// storage of intermediate results
		struct TensorInfo[2] tensmem; // tensor memory
		// backup tensor memory for branching tensor calculations
		struct TensorInfo[2] tensbmem;
		// indices of the last optimization center used
		// the indices are stored such that last_optimized
		int last_optimized[2];
};


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
	// prepare the symsecs of the result
	result->renew_symsec_layout(reference, OO_link);
	// prepare the block layout of the result
	result->renew_block_layout(reference, true);

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
	// assert that OO is associated to the open leg
	assert(OO->ref == ref->syms[open_leg] && OO->opt == opt->syms[open_leg])

	// Indices are permuted (=:per) such that the outermost loop corresponds
	// with the index in the open leg == block index in the OverlapObject
	int per[3] = {open_leg, (open_leg+1)%3,(open_leg+2)%3}
	// Collect the matching irrep pairs for each leg 
	for (int i=0; i<3; i++) {
		symsecMatchers[i].set_matching_symsec_indices(&ref, &opt, i);
	}

	// modify the block layout of OO
	OO->renew_block_layout(&symsecMatchers[leg] , true)

	// irrep indices for the reference and optimizing blocks currently
	// contributing
	int r[3], o[3];

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

/*****************************SymsecMatcher****************************/

class SymsecMatcher {
	public:
		// constructor and destructor
		SymsecMatcher(int maxSize) : size(0) {
			result = safe_malloc(2 * maxSize, int); }
		~ SymsecMatcher() { free(result); }

		// get the current results
		int get_size() const { return size; };
		const int ** get_result() const { return result; };

		// calculate the new results
		// The results are guaranteed to be sorted in the irrep index
		// of the first argument.
		void set_matching_symsec_indices(const TensorInfo * a,
			const TensorInfo * b, int leg) const;
		void set_matching_symsec_indices(const struct symsecs * a,
			const struct symsecs * b) const;
	private:
		// data fields
		int size;  // effective space currently used
		int (*result)[2];  // pointer to array of int pairs
};

void SymsecMatcher::set_matching_symsec_indices(const TensorInfo * A,
			const TensorInfo * B, int leg) const {
	get_matching_symsec_indices(A->symsecs[leg], B->symsec[leg]);
}

void SymsecMatcher::set_matching_symsec_indices(const struct symsecs * a,
			const struct symsecs * b) const {
	int index = 0;
	for (int i = 0; i < a->nrSecs; i++) {
		j = search_symsec(&a.irreps[i], b);
		if (j != -1) {
			(this->result)[index] = {i,j};
			index++;
		}
	}
	this->size = index;
}

/*********************************TensorInfoPair**************************/

class TensorInfoPair {
	// tensor of the reference state, e.g. the ground state
	TensorInfo ref;
	// tensor of the state that is being optimized,
	// e.g. the excited state
	TensorInfo opt;
};

class TensorInfo {
	public:
		// constructor
		TensorInfo(siteTensor * data, struct symsecs ** syms);

		// Get the index of the OO-block corresponding to the block
		// with given index of the data when the OO is associated
		// with the bond 'leg'.
		int get_leg_index(int index, int leg);
		void get_leg_indices(int index, int * result);

		// Get the contraction info for the block with given indices
		// ->  labeling with irrep indices of the bonds
		struct BlockInfo get_block_info(int i, int j, int k);
		// ->  labeling with the indices from in data
		struct BlockInfo get_block_info(int index);

		// renew the symmetry sectors based
		// @param:
		//   ref => contains the uncontracted symsecs
		//   OO_link => maps the contracted symsec onto the optimizing one
		void renew_symsec_layout(const TensorInfo * ref,
			const OverlapObjectLink * OO_link)

		// renew all block information; optionally set elements to 0
		// @param:
		//   reference => contains the qnumbers
		//   set_zero = fill the tel array with zeros
		void renew_block_layout(const TensorInfo * reference, bool set_zero)

		// set result equal to the contraction of data with toAdd
		void add_up_with_OverlapObject(OverlapObjectLink * toAdd,
				TensorInfo * result, TensorInfo * memory);

		// get or set the symmetry sectors
		struct symsecs ** get_syms() const { return syms; }
		struct symsecs * get_sym(int leg) const { return syms[leg]; }
		void set_syms(struct symsecs ** new_syms) {
			for (int i=0; i<3; i++) { syms[i] = new_syms[i]; }}
		void set_sym(struct symsecs * new_sym, int leg) {
			syms[leg] = new_sym; }

		// is the tensor a physical or branching tensor?
		bool is_physical_tensor() const { return is_physical; }
		// get the data
		struct siteTensor * get_data() const { return data; }
	private:
		// symmetry structures in order |a>|b><c| (cf. siteTensor)
		struct symsecs * syms[3];
		// actual siteTensor containing the data
		struct siteTensor * data;
		// is the tensor a physical or branching tensor?
		bool is_physical;
		// size of tel array used to effectively store elements
		int usedsize;
		// allocated size of the tel array; should be >= usedsize
		int allocsize;
};


void renew_symsec_layout(const TensorInfo * ref,
			const OverlapObjectLink * OO_link)
{
	// extract the index of the leg over which is contracted
	int contracted_leg = OO_link->leg;

	// check whether the contraction is valid
	assert(OO_link->OO->get_ref() == reference->get_sym(OO->leg));

	// set the relevant sectors
	result->set_syms(ref->get_syms());
	result->set_sym(contracted_leg, OO_link->OO->get_opt());
}


void renew_block_layout(const TensorInfo * ref, bool set_zero)
{
	// number of blocks of the new structure
	int nrblocks = ref->get_data()->nrblocks;

	// reallocate qnumbers and beginblock if necessary
	if (data->nrblocks < nrblocks) {
		realloc(data->qnumbers, nrblocks, sizeof(QN_TYPE));
		realloc(data->blocks->beginblock, nrblocks, sizeof(int));
		nr_allocated_blocks = nrblocks;
	}

	// set number of blocks
	data->nrblocks = nrblocks;

	// fill qnumbers
	for (int i=0; i<data->nrblocks; i++) {
		data->qnumbers[i] = ref->data->qnumbers[i];
	}

	// fill beginblock
	int usedsize = 0;
	int dims[3];  // dimensions of the current block
	for (int i=0; i<data->nrblocks; i++) {
		get_block_dimensions(i, dims);
		data->blocks->beginblock[i] = usedsize;
		usedsize += dims[0] * dims[1] * dims[2];
	}

	// reallocate tel if necessary
	if (usedsize > allocsize) {
		allocsize = 2 * usedsize;
		this->reallocate(allocsize); }
	if (set_zero) {
		// fill tel up with zeros
		for (int i=0; i<usedsize; i++) {
			blocks.tel[i] = 0; }
	}
}


inline int get_block_dimension(int leg_index){
	return (syms[i])->dims[leg_index]; }


inline void TensorInfo::get_block_dimensions(int index, int * dims)
{
	// get the leg indices
	int ids[3];
	get_leg_indices(index, ids);
	// extract the corresponding dimensions
	for (int i=0; i<3; i++) {
		dims[i] = get_block_dimension(ids[i]); }
}


inline int TensorInfo::get_leg_index(int index, int leg)
{
	int ids[3];
	get_leg_indices(index, ids);
	return ids[leg];
}


inline void TensorInfo::get_leg_indices(int index, int * result) {
	indexize(result, data->qnumbers[index], syms); }


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


// // WARNING: WRONG PREFACTORS & ALMOST TRIVIAL ARE POSSIBLE BY USING GAUGE FREEDOM
// double get_2leg_contraction_prefactor(int open_leg, TensorInfo reference,
// 	int * irrep_indices)
// {
// 	// irrep indices are renamed such that o refer to the open leg,
// 	// a to the 1st leg following o and b to the 2nd leg following o
// 	int o = open_leg;        // index of the open leg
// 	int a = (open_leg+1)%2;  // first leg following o
// 	int b = (open_leg+2)%2;  // second leg following o
// 	// similar naming conventions are applied for the total angular
// 	// momentum (j) and parity (p) associated with them
// 	int jo, ja, jb, po, pa, pb;

// 	// hold a pointer to the relevant symsecs
// 	struct symsecs * syms = reference->get_syms();

// 	jo = get_angular_momentum(irrep_indices[o], syms[o]);
// 	jb = get_angular_momentum(irrep_indices[b], syms[b]);
// 	pa = get_parity(irrep_indices[a], syms[a]);

// 	// ALARM: check correctness
// 	if (reference->is_physical_tensor()) {
// 		assert(open_leg == 0);
// 		result = grade(pp) * (2 * jb + 1) / (2 * jo + 1);
// 	} else {
// 		ja = get_angular_momentum(irrep_indices[a], syms[a]);
// 		pb = get_parity(irrep_indices[b], syms[b]);
// 		switch (open_leg)
// 		{
// 			case 0, 1:
// 				po = get_parity(irrep_indices[o], syms[o]);
// 				jx = (open_leg == 0)? ja : jb;
// 				return grade(ja+jb+jo + pa+pb+po) *
// 					(2 * jx + 1) / (2 * jo + 1);
// 					// ALARM: pa = pb*po?
// 				break;

// 			case 2:
// 				result = grade(pa*pb + ja+jb+jo);
// 				break;
// 		}
// 	}
// }


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


/*********************************OverlapObject****************************/

// IDEA: DO NOT REMEMBER REF AND OPT BUT TAKE IT AS INPUT

// Usage: ones created, first renew the block layout with a solved
// symsecMatcher before performing block calculations that take block info.
class OverlapObject {
	public:
		// constructor
		// @param:
		//   ref = symmetry sectors of the bond in the reference T3NS
		//   opt = symmetry sectors of the bond in the optimizing T3NS
		//   [opt_dim] = current/bond dimension of opt
		// @remark: if opt_dim is set to the bond dimension, no reallocation
		//   will be needed. Remark that 2000**2 * 100 * 8 bytes is about 3Gb.
		OverlapObject(struct symsecs * ref, struct symsecs * opt,
				int opt_dim = 2*opt->totaldims);

		// destructor
		~ OverlapObject() { free(blocks.beginblock); free(blocks.tel);
				free(ldim); free(sdim); }

		// renew all block information; optionally set elements to 0
		// @param:
		//   match => contains the matching irrep indices between ref and opt
		//   set_zero = fill the tel array with zeros
		void renew_block_layout(const SymsecMatcher * match, bool set_zero)

		// get the info for a block contraction
		// Remark that result.start == -1 if the block doesn't exists.
		struct BlockInfo get_block_info(int index) const;

		// reallocate the element array explicitely
		// @param: opt_dim = total/bond dimension for the optimizing bond
		void reallocate(int opt_dim);

		// get the reference or optimizing symsec
		struct symsec * get_ref() { return ref; }
		struct symsec * get_opt() { return opt; }

	private:
		// Symmetry sectors in the reference and optimizing bond
		// to wich the OO is associated
		struct symsecs * ref, * opt;
		// Actual elements of the OverlapObject.
		//   The block corresponding to the i'th irrep of ref starts at
		//   blocks.tel[blocks.beginblock[i]], has leading dimension
		//   ldim[i] and second dimension sdim[i].
		//   If blocks.beginblock[i] == -1, the i'th block is not present.
		struct sparseblocks blocks;  // begin blocks and element array
		int * ldim;  // leading dimension array
		int * sdim;  // second dimension array
		// size of tel array used to effectively store elements
		int usedsize;
		// allocated size of the tel array; should be >= usedsize
		int allocsize;
};


OverlapObject::OverlapObject(struct symsecs * ref, struct symsecs * opt,
			int opt_dim = 2*opt->totaldims) : ref(ref), opt(opt),
			telsize(opt_dim)
{
	// help variable
	int nrBlocks = ref->nrSecs;

	// initialize the sparseblocks
	init_null_sparseblocks(blocks);
	// -> allocate the tel array
	blocks.tel = safe_malloc(ref->totaldims * opt_dim, EL_TYPE);
	// -> allocate beginblock
	blocks.beginblock = safe_malloc(nrBlocks, int);

	// allocate the ldim and sdim array
	ldim = safe_malloc(nrBlocks, int);
	sdim = safe_malloc(nrBlocks, int);
}


// renew the beginblocks and tel array
// match should contain the matches between ref and opt
// no assumptations are made about the ordering of matching symsecs
void OverlapObjects::renew_block_layout(const SymsecMatcher * match,
		bool set_zero)
{
	// help variables
	usedsize = 0;
	int nrBlocks = ref->nrSecs;
	// reset the beginblock and dimension arrays
	for (int i=0; i<nrBlocks+1; i++) {
		blocks.beginblock[i] = -1;  // default value
		ldim[0] = sdim[0] = 0;
	}
	// fill the beginblock and dimension arrays
	for (int j=0; j<match.get_size(); j++) {
		// get the relevant indices
		int ref_index = match.get_result()[j][0];
		int opt_index = match.get_result()[j][1];
		// set the begin and dims of the block
		beginblock[ref_index] = usedsize;
		ldim[ref_index] = ref->dims[ref_index];
		sdim[opt_index] = opt->dims[opt_index];
		// increase the usedsize
		usedsize += ldim[ref_index] * sdim[opt_index];
	}
	// reallocate tel if necessary
	if (usedsize > allocsize) {
		allocsize = 2 * usedsize;
		this->reallocate(allocsize); }
	if (set_zero) {
		// fill tel up with zeros
		for (int i=0; i<usedsize; i++) {
			blocks.tel[i] = 0; }
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


// reallocate when the used dimension for the optimizing T3NS changes
void OverlapObject::reallocate(int opt_dim)
{
	blocks.tel = realloc(blocks.tel, ref->totaldims * max_opt_dim *
			sizeof(EL_TYPE));
}


struct OverlapObjectLink {
	OverlapObject* OO;
	int leg;  // index of the leg to which the OOlink is associated
};





// TRASH


// class PhysicalTensorInfo {
// 	public:
// 		void add_up_with_OverlapObject(const OverlapObjectLink * toAdd,
// 			TensorInfo * result, TensorInfo * memory);
// };

// void PhysicalTensorInfo::add_up_with_OverlapObject(const OverlapObjectLink * toAdd,
// 		TensorInfo * result, TensorInfo * memory) {
// 	for (int i=0; i < data->nrBlocks; i++) {
// 		double prefactor = getPrefactor(); //TODO
// 		j = get_OO_index(i,toAdd->leg);
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [result]_ijk = prefactor * [this]_ijl * [OO]_lk
// 		block_product(prefactor,
// 			this->get_block_info(i),
// 			toadd->OO->get_block_info(j),
// 			result->get_block_info(i));
// 	}
// }

// class BranchingTensorInfo {
// 	public:
// 		void add_up_with_OverlapObject(const OverlapObjectLink * toAdd,
// 			TensorInfo * result, TensorInfo * memory);
// }

// void BranchingTensorInfo::add_up_with_OverlapObject(const OverlapObjectLink * toAdd,
// 		TensorInfo * result, TensorInfo * memory) {
// 	for (int i=0; i < data->nrBlocks; i++) {
// 		double prefactor = getPrefactor(); //TODO
// 		j = get_OO_index(i,toAdd[0]->leg);
// 		k = get_OO_index(i,toAdd[1]->leg);
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [memory]_ijk = prefactor * [this]_ijl * [OO[0]]_lk
// 		// -> memory = this * OO[0]
// 		block_product(1,
// 			this->get_block_info(i),
// 			toadd[0]->OO->get_block_info(j),
// 			memory->get_block_info(i));
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [result]_ijk = prefactor * [memory]_ijl * [OO[1]]_lk
// 		// -> result = prefactor * memory * OO[1]
// 		block_product(prefactor,
// 			memory->get_block_info(i),
// 			toadd[1]->OO->get_block_info(k),
// 			result->get_block_info(i));
// 	}
// }

// oid set_OO_block_to_contraction(double prefactor, int leg,
// 		const struct BlockInfo A, const struct BlockInfo B,
// 		struct BlockInfo C) {
// 	struct contractinfo cinfo = {
// 		.tensneeded = {0, 1, 2},
// 		. = {CblasNoTrans, CblasNoTrans},
// 		// TODO
// 		.M = A.d1,
// 		.N = B.d3,
// 		.K = A.,
// 		.L = ,
// 		.lda = ,
// 		.ldb = ,
// 		.ldc = ,
// 		.stride = 
// 	};
// 	EL_TYPE tels[3] = {A->start, B->start, C->start};
// 	if (leg == 0) {
// 		cblas_dgemm()
// 	}
// }

// // WARNING: wrong formula: C_ij = prefactor * A_ij * B_jk
// void block_product(double prefactor, const struct BlockInfo A,
// 		const struct BlockInfo B, struct BlockInfo C) {
// 	struct contractinfo cinfo = {
// 		.tensneeded = {0, 1, 2},
// 		. = {CblasNoTrans, CblasNoTrans},
// 		// TODO
// 		.M = A.,
// 		.N = B.,
// 		.K = A.,
// 		.L = ,
// 		.lda = ,
// 		.ldb = ,
// 		.ldc = ,
// 		.stride = 
// 	};
// 	EL_TYPE tels[3] = {A->start, B->start, C->start};
// 	do_contract(cinfo, tels, prefactor, 1);
// }

// // perform cblas_dgemm, but implicitely pad zeros to match dimensions
// void pad_zero_dgemm(CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
// 		int M, int N, int K, double alpha, double * A, int lda,
// 		double * B, int ldb, double beta, double * C, int ldc)
// {

// }

//void add_Aikl_C_ij_eq

// blocks.beginblock[0] = counter;
// 	match.get_result()[i][0]
// 	for (int i=0; i<nrBlocks; i++) {
// 		counter += ref->dims[i] * opt->dims[i];
// 		blocks.beginblock[i+1] = counter; }
// 	// reallocate tel if necessary
// 	if (counter > telsize) {
// 		telsize = 2 * counter;
// 		this->reallocate(telsize); }
// 	if (set_zero) {
// 		// fill tel up with zeros
// 		int length = blocks.beginblock[nrBlocks];
// 		for (int i=0; i<length; i++) {
// 			blocks.tel[i] = 0; }
// 	}

// // tensor A is considered as the reference for OO's symsectioning
// void TwoSiteOverlapCalculator::set_OO_to_contraction(const TensorInfo * A,
// 			const TensorInfo * B, int leg, OverlapObject * OO) {
// 	OO->reset_blocks();
// 	for (int i=0; i<3; i++) {
// 		symsecMatchers[i].set_matching_symsec_indices(&A, &B, i);
// 	}
// 	int a[3], b[3];
// 	for (int i=0; i<symsecMatchers[0].get_size(); i++) {
// 		a[0] = symsecMatchers[0].get_result()[i][0];
// 		b[0] = symsecMatchers[0].get_result()[i][1];
// 		for (int j=0; j<symsecMatchers[1].get_size(); i++) {
// 			a[1] = symsecMatchers[1].get_result()[j][0];
// 			b[1] = symsecMatchers[1].get_result()[j][1];
// 			for (int k=0; k<symsecMatchers[2].get_size(); i++) {
// 				a[2] = symsecMatchers[2].get_result()[k][0];
// 				b[2] = symsecMatchers[2].get_result()[k][1];
// 				// do the actual block manipulation
// 				// ALARM: ORDER OF MULTIPLICATION???
// 				// [OO]_ij = prefactor * A_kli * B_klj

// 				// if with RG-flow
// 				// [OO]_block(k) = A_block(i,j,k) * B_block(i,j,k)

// 				// if against RG-flow
// 				// [OO]_block(i) = A_block(i,j,k) * B_block(i,j,k)
// 				block_product(get_prefactor(),
// 					A->get_block_info(a[0],a[1],a[2]),
// 					B->get_block_info(b[0],b[1],b[2]),
// 					OO->get_block_info(a[leg]));
// 			}
// 		}
// 	}
// }


// OverlapObject::OverlapObject(struct * symsecs ref, struct * symsecs opt,
// 			bool greedy) : ref(ref), opt(opt), greedy(greedy) {
// 	// help variables
// 	int counter, nrBlocks, opt_dim;
// 	nrBlocks = ref->nrSecs;

// 	// initialize the sparseblocks
// 	init_null_sparseblocks(blocks);
// 	opt_dim = (greedy)? max_opt_dim : opt->totaldims;
	
// 	blocks.tel = safe_malloc(ref->totaldims * opt_dim, EL_TYPE);


// 	// -> construct the beginblock array
	
	
// 	blocks.beginblock = safe_malloc(nrBlocks, int);
// 	blocks.beginblock[0] = counter = 0;
// 	for (int i=0; i<nrSecs; i++) {
// 		counter += ref->dims[i] * opt->dims[i];
// 		blocks.beginblock[i+1] = counter; }
// 	// -> allocate the tel array
// 	blocks.tel = safe_malloc(counter, EL_TYPE);
// }

// // If greedy, (D_ref x D_opt x 8) bytes was allocated and no
		// // reallocation is needed during recalculation.
		// // Remark that 2000**2 * 100 * 8 bytes is about 3Gb.
		// bool greedy;

// void OverlapObject::reset_blocks() {
// 	int length = blocks.beginblock[syms.nrSecs];
// 	for (int i=0; i<length; i++) {
// 		blocks.tel[i] = 0;
// 	}
// }

// do the actual block manipulation
				// ALARM: ORDER OF MULTIPLICATION???
				// [OO]_ij = prefactor * A_kli * B_klj

				// if 'leg = 0' (optimization centre moves against RG-flow):
				// [OO]_block(i) += prefactor * A_block(i,j,k) * B_block(i,j,k)
				// if 'leg = 2' (optimization centre moves with RG-flow):
				// [OO]_block(k) += prefactor * A_block(i,j,k) * B_block(i,j,k)

// for (int j=0; j<2; j++) {
// 			if (optimization_center[i] == last_optimized[j]) {
// 				common = 
// 			}
// 		}
// 	}
// 	assert (linsearch(opt_center[0], last_optimized, 2, sort_int[1], sizeof(int)) != -1
// 		|| linsearch(opt_center[0], last_optimized, 2, sort_int[1], sizeof(int)) != -1);


// TwoSiteOverlapCalculator::get_overlap_vector(int * optimization_center,
// 			double * result) {
// 	for (int i=0; i<2; i++) {
// 		int f = optimization_center[i]    // first tensor index
// 		int l = optimization_center[1-i]  // last tensor index
// 		TensorInfoPair current = tensorpairs[f];
// 		// do not copy the content of the data blocks
// 		shallow_copy(current->ref, tensmem[i]);
// 		OverlapObjectLink external[2] = external_links[f,l];
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [current->ref]_ijk = [external->OO]_il * [tensmem]_ljk
// 		current->ref->add_up_with_OverlapObject(external,
// 				tensmem[i], tensbmem[i]);
// 		// TODO
// 		// "change_sector_tensors"(tensmem[i], current->opt, all external->leg)
// 	    // => match the tensmem[i]->symsec[external->leg]
// 	    //    with current->opt->symsec[external->leg] for all legs
// 	    // => change tensmem[i] to a vector vec_x
// 	    // remember the calculation of tensmem[i]
// 	    last_optimized[i] = f;
// 	}
// 	// TODO
// 	makeSiteTensor(vec_a, vec_b);
// }

// The tensor for which the OO needs to be updates was also part of the previous
	// optimization center. The other tensor indicates the bond which is NOT involved
	// in update of the


// // The optimizing tensors are assumed to correspond to the indices of the network
// // used in the constructor. No further assumptions are made.
// TwoSiteOverlapCalculator::get_overlap_vector(int * optimizing_tensors, double * result)
// {
// 	// First, the Overlap Object between the previously updated sites might be updated
// 	// so that it can be used in the calculation of the new overlap vector. This means
// 	// that the particular OO needs to be directed to the current optimization center.
// 	int * common = null;      // common index with previous optimization center
// 	int * left_over = null;   // index that's only in the previous optimization center
// 	TensorInfo * memory;      // memory location of the intermediate saved result

// 	// compare the current and previous optimization center
// 	for (int i=0; i<2; i++) {
// 		int& current = last_optimized[i];
// 		if(linSearch(current, optimizing_tensors, 2, sort_int[1], sizeof(int)) == -1) {
// 			memory = &(tensmem[i]);
// 			*left_over = current;
// 		} else {
// 			*common = current; }
// 	}

// 	// A necessary condition is that the current and previous optimization
// 	// center have at least one tensor in common.
// 	assert(common != null);

// 	// Update if the optimization center has moved
// 	if (left_over != null) {
// 		update(*left_over, *common, memory);
// 	}

// 	// Second, ....
// 	// for each of the tensors involved
// 	for (int i=0; i<2; i++) {
// 		int f = optimizing_tensors[i]    // first tensor index
// 		int l = optimizing_tensors[1-i]  // last tensor index
// 		TensorInfoPair current = tensorpairs[f];
// 		// do not copy the content of the data blocks
// 		shallow_copy(current->ref, tensmem[i]);
// 		OverlapObjectLink external[2] = external_links[f,l];
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [current->ref]_ijk = [external->OO]_il * [tensmem]_ljk
// 		current->ref->add_up_with_OverlapObject(external,
// 				tensmem[i], tensbmem[i]);
// 		// TODO
// 		// "change_sector_tensors"(tensmem[i], current->opt, all external->leg)
// 	    // => match the tensmem[i]->symsec[external->leg]
// 	    //    with current->opt->symsec[external->leg] for all legs
// 	    // => change tensmem[i] to a vector vec_x
// 	    // remember the calculation of tensmem[i]
// 	    last_optimized[i] = f;
// 	}
// 	// TODO
// 	makeSiteTensor(vec_a, vec_b);

// }

// TwoSiteOverlapCalculator::update() {
// 	if (last_optimized[0] == -1) {
// 		std::cout << "ALARM!";
// 	} else {
// 		int f = last_optimized[0];  // first tensor index
// 		int l = last_optimized[1];  // last tensor index
// 		TensorInfoPair current = tensorpairs[f];
// 		OverlapObjectLink * internal = internal_link[f,l];
// 		// Save in internal->OO the contraction of tensmem and the relevant
// 		// site of the optimizing T3NS.
// 		set_OO_to_contraction(&tensmem[0], &(current->opt),
// 				internal->leg, internal->OO);
// 		last_optimized[0] = -1;
// 	}
// }

// // START TO DO REGION
// 		// do not copy the content of the data blocks
// 		shallow_copy(current_pair->ref, tensmem[i]);
// 		OverlapObjectLink external[2] = external_links[f,l];
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [current->ref]_ijk = [external->OO]_il * [tensmem]_ljk
// 		current->ref->add_up_with_OverlapObject(external,
// 				&(tensmem[i]), &(tensbmem[i]));
// 		// END TO DO REGION


// void OverlapCalculator::contract_tensor_with_OO_list(const TensorInfo * tensor,
// 			const OverlapObjectLink * OO_link_list, TensorInfo * memory,
// 			TensorInfo * result)
// {
// 	if (tensor->is_physical_tensor()) {
// 		for (int i=0; i < data->nrBlocks; i++) {
// 			double prefactor = getPrefactor(); //TODO
// 			j = get_OO_index(i,toAdd->leg);
// 			// ALARM: ORDER OF MULTIPLICATION???
// 			// [result]_ijk = prefactor * [this]_ijl * [OO]_lk
// 			block_product(prefactor,
// 				this->get_block_info(i),
// 				toadd->OO->get_block_info(j),
// 				result->get_block_info(i));
// 		}
// 	} else {
// 	for (int i=0; i < data->nrBlocks; i++) {
// 		double prefactor = getPrefactor(); //TODO
// 		j = get_OO_index(i,toAdd[0]->leg);
// 		k = get_OO_index(i,toAdd[1]->leg);
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [memory]_ijk = prefactor * [this]_ijl * [OO[0]]_lk
// 		// -> memory = this * OO[0]
// 		block_product(1,
// 			this->get_block_info(i),
// 			toadd[0]->OO->get_block_info(j),
// 			memory->get_block_info(i));
// 		// ALARM: ORDER OF MULTIPLICATION???
// 		// [result]_ijk = prefactor * [memory]_ijl * [OO[1]]_lk
// 		// -> result = prefactor * memory * OO[1]
// 		block_product(prefactor,
// 			memory->get_block_info(i),
// 			toadd[1]->OO->get_block_info(k),
// 			result->get_block_info(i));
// 	}
// }
// }


// OverlapObjectLink to_prepare[2];
// 		nr_to_prepare = get_external_links(who, exclude, to_prepare);
// 		for (int i=0; i<nr_to_prepare; i++) {
// 			if (to_prepare[i] == -1) {
// 				get_internal_link(who, exclude,)
// 			}
// 		}
// 	}


// for (int i=0; i<size; i++) {

// 		}


// 		// ... look up their adresses, ...
// 		get_bonds_of_site(who, bond_inds);
// 		// ... ask them to do their part ...
// 		// - get the bond indices for 'who'
// 		get_
// 		// - loop over the bond indices
// 		for (int i=0; i<3; i++) {
// 			// - if the bond is a virtual bond
// 			if (bond_inds[i] != -1) {
// 				// - get the pairs linked by the bond
// 				bond_link = &(network->bonds[bond_inds[i]]);
// 				int from = (*bond_link)[0];
// 				int to = (*bond_link)[1];
// 				// - if not excluded, prepare the OOs for the given neighbour
// 				if (from != exclude && to != exclude) {
// 					prepare_OOs_for_pair(who, to); }}}

// 		// ... and merge their results.
// 		for (int i=0)


// int OverlapCalculator::get_external_links(int who, int exclude,
// 		OverlapObjectLink * result)
// {
// 	int bond_inds[3];     // bond indices
// 	int (*bond_link)[2];  // adres of current bond
// 	int size = 0;         // size of the result

// 	// get the bond indices for 'who'
// 	get_bonds_of_site(who, bond_inds);
// 	// loop over the bond indices
// 	for (int i=0; i<3; i++) {
// 		// if the bond is a virtual bond
// 		if (bond_inds[i] != -1) {
// 			// ask the tensors linked by the bond
// 			bond_link = &(network->bonds[bond_inds[i]]);
// 			int from = (*bond_link)[0];
// 			int to = (*bond_link)[1];
// 			// if not excluded, add an OO_link to the result
// 			if (from != exclude && to != exclude) {
// 				result[size].OO = &overlaps[bond_inds[i]];
// 				result[size].leg = i;
// 				size++; }}}

// 	// return the size of the final result
// 	return size;
// }


// // the indexing of result corresponds to that of site tensors in opt
// void OverlapCalculator::get_neighbouring_pairs(int who, int * result)
// {
// 	int bond_inds[3];     // bond indices
// 	int (*bond_link)[2];  // adres of current bond

// 	// get the bond indices for 'who'
// 	get_bonds_of_site(who, bond_inds);
// 	// loop over the bond indices
// 	for (int i=0; i<3; i++) {
// 		// if the bond is a virtual bond
// 		if (bond_inds[i] != -1) {
// 			// get the tensors linked by the bond
// 			bond_link = &(network->bonds[bond_inds[i]]);
// 			int from = (*bond_link)[0];
// 			int to = (*bond_link)[1];
// 			// append the other pair index to the result
// 			result[i] = (from == who)? to : from;
// 		} else {
// 			result[i] = -1; }}
// }


// get_other_pair


// void OverlapCalculator::get_internal_link(int who, int other,
// 		OverlapObject * result)
// {
// 	// collect the indices of the neighbouring pairs
// 	int neighbours[3];
// 	get_neighbouring_pairs(who, neighbours);
// 	for (int i=0; i<size; i++) {
// 		if (neighbours[i] != -1) {
// 			result->OO = &(overlaps[bond_inds[i]]);
// 		}
// 		// if not excluded, add an OO_link to the result
// 			if (to == other) {
// 				result->OO = &overlaps[bond_inds[i]];
// 				result->leg = i; }}}
// 	}
// }


// void OverlapCalculator::get_internal_link(int who, int other,
// 		OverlapObject * result)
// {
// 	int bond_inds[3];     // bond indices
// 	int (*bond_link)[2];  // adres of current bond
// 	int size = 0;         // size of the result

// 	// get the bond indices for 'who'
// 	get_bonds_of_site(who, bond_inds);
// 	// loop over the bond indices
// 	for (int i=0; i<3; i++) {
// 		// if the bond is a virtual bond
// 		if (bond_inds[i] != -1) {
// 			// ask the tensors linked by the bond
// 			bond_link = &(network->bonds[bond_inds[i]]);
// 			int from = (*bond_link)[0];
// 			int to = (*bond_link)[1];
// 			// if not excluded, add an OO_link to the result
// 			if (to == other) {
// 				result->OO = &overlaps[bond_inds[i]];
// 				result->leg = i; }}}
// }


// bool with_physical_tensor,
// 	int open_leg, int * irrep_indices, const struct symsecs * syms)

// double get_2leg_contraction_prefactor(int open_leg, TensorInfo reference,
// 	int * irrep_indices)
// {
// 	// irrep indices are renamed depending on the coupling in reference
// 	int o = irrep_indices[open_leg];  // irrep index of the open leg
// 	int p;  // irrep index of the physical leg
// 	int a;  // irrep index of the first branching leg in the coupling
// 	int b;  // irrep index of the second branching leg in the coupling
// 	// similar naming conventions are applied for the total angular
// 	// momentum (j) and parity (p) associated with it
// 	int ja, jb, jo, pp, pa, pb, po;

// 	// hold a pointer to the relevant symsecs
// 	struct symsecs * syms = reference->get_syms();

// 	if (reference->is_physical_tensor()) {
// 		assert(open_leg == 0);
// 		p = 1; a = 2;
// 		// get pp, ja and jo
// 		pp = get_parity(irrep_indices[p], syms[p]);
// 		ja = get_angular_momentum(irrep_indices[a], syms[a]);
// 		jo = get_angular_momentum(irrep_indices[o], syms[o]);
// 	} else {
// 		switch (open_leg)
// 		{
// 			case 0:
// 				a = irrep_indices[1];
// 				b = irrep_indices[2];
// 				// get pp, ja and jo
// 				break;

// 			case 1, 2:
// 				a = irrep_indices[0];
// 				b = irrep_indices[3 - open_leg];
// 				// get ja, jb, jo, pb and po
// 				break;
// 		}
// 	}
// 	// do magic
// }

// double get_prefactor() {
// 	// TODO
// 	std::cout << "this function should be implemented correctly" << endl;
// 	std::cout << "the right parameters should be asked" << endl;
// 	return 1;
// }

// if (opt_t3ns->nr_data != ref_t3ns->nr_data) {
// 		throw invalid_argument("Number of sites in the reference " +
// 				"T3NS and optimizing T3NS should match."); }
// 	if (opt_t3ns->nr_data != netw->sites) {
// 		throw invalid_argument("Number of sites in the network and the " +
// 				"optimizing T3NS should match."); }
// 	if (opt_t3ns->bookie->nr_bonds != ref_t3ns->bookie->nr_bonds) {
// 		throw invalid_argument("Number of bonds in the reference " +
// 				"T3NS and optimizing T3NS should match."); }
// 	if (opt_t3ns->bookie->nr_bonds != netw->nr_bonds) {
// 		throw invalid_argument("Number of bonds in the network and the " +
// 				"optimizing T3NS should match."); }

// returns 0 (even parity) or 1 (odd parity)
// int get_parity(int irrep_index, const struct * symsecs syms)
// {
// 	char buffer[MY_STRING_LEN];
// 	int result;
// 	get_irrstring(buffer, Z2, syms->irreps[irrep_index]);
// 	which_irrep(buffer, Z2, &result);
// 	return result;
// }


// // returns 2j+1?
// // ALARM: IS THIS REALLY THE WAY TO DO IT?
// double get_angular_momentum(int irrep_index, const struct * symsecs syms)
// {
// 	char buffer[MY_STRING_LEN];
// 	int result;
// 	get_irrstring(buffer, SU2, syms->irreps[irrep_index]);
// 	which_irrep(buffer, SU2, &result);
// 	return result;
// }

// // returns 0 (even parity) or 1 (odd parity)
// int get_angular_momentum(int irrep_index, const struct symsecs * syms,
// 		const struct bookkeeper * bookie)
// {
// 	for (int i = 0; i < bookie.nrSyms; i++) {
// 	    if (bookie.sgs[i] == Z2) {
// 	        return syms->irreps[irrep_index][i]; }}

// 	// No Z2 is present.
// 	fprintf(stderr, "Warning: no Z2 is used.\n");
// 	return -1;
// }


// // returns 2j
// int get_angular_momentum(int irrep_index, const struct symsecs * syms,
// 		const struct bookkeeper * bookie)
// {
// 	for (int i = 0; i < bookie.nrSyms; i++) {
// 	    if (bookie.sgs[i] == SU2) {
// 	        return syms->irreps[irrep_index][i]; }}

// 	// No SU2 is present.
// 	fprintf(stderr, "Warning: no SU2 is used.\n");
// 	return -1;
// }

#ifndef TENSOR_INFO
#define TENSOR_INFO

#include "siteTensor.h"
#include "symsecs.h"


class TensorInfo {
	public:
		// constructor
		TensorInfo(struct siteTensor * data, struct symsecs ** syms) {};
	private:
		// Main data:
		// symmetry structures in order |a>|b><c| (cf. siteTensor)
		struct symsecs (*syms)[3];
		// actual siteTensor containing the data
		struct siteTensor * data;

		// Help data:
		// is the tensor a physical or branching tensor?
		bool is_physical;
		// size of tel array used to effectively store elements
		int usedsize;
		// allocated size of the tel array; should be >= usedsize
		int allocsize;
};


typedef struct TensorInfoPair {
	// tensor of the reference state, e.g. the ground state
	TensorInfo ref;
	// tensor of the state that is being optimized,
	// e.g. the excited state
	TensorInfo opt;
} TensorInfoPair;


#endif


#ifndef TENSOR_INFO
#define TENSOR_INFO

#include "cpp_interface.h"
#include "siteTensor.h"
#include "symsecs.h"
#include <stdio.h>
#include <iostream>


class TensorInfo {
	public:
		// constructor
		TensorInfo(struct siteTensor * data, struct symsecs ** syms);

		// get or set the symmetry sectors
	//	struct symsecs ** get_syms() const { return &syms[0]; }
		struct symsecs * get_sym(int leg) const { return syms[leg]; }
		void set_syms(struct symsecs ** new_syms) {
			for (int i=0; i<3; i++) { syms[i] = new_syms[i]; }}
		void set_sym(struct symsecs * new_sym, int leg) {
			syms[leg] = new_sym; }

		// get the data
		struct siteTensor * get_data() const { return data; }
	private:
		// Main data:
		// symmetry structures in order |a>|b><c| (cf. siteTensor)
		struct symsecs * syms[3];
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


void print_tensorinfo(const TensorInfo * tensor,
		const struct bookkeeper * keeper);
void print_tensorinfopair(const TensorInfoPair * pair,
		const struct bookkeeper * keeper);
// std::ostream& operator << (std::ostream& out, const TensorInfo& tensor);


#endif

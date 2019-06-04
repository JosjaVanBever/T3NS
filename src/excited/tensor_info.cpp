
#include "tensor_info.h"

TensorInfo::TensorInfo(struct siteTensor * data, struct symsecs ** syms) :
		data(data)
{
	for (int i=0; i<3; i++) { this->syms[i] = syms[i]; }
}


void print_tensorinfo(const TensorInfo * tensor, const struct bookkeeper * keeper)
{
	struct siteTensor * data = tensor->get_data();
	print_siteTensor(keeper, data);
	// struct symsecs * syms[3] = tensor->get_sym(i);
	for (int i=0; i<3; i++) {
		struct symsecs * sym = tensor->get_sym(i);
		print_symsecinfo(sym);
	}
}

void print_tensorinfopair(const TensorInfoPair * pair,
		const struct bookkeeper * keeper)
{
	print_tensorinfo(&(pair->ref), keeper);  // WRONG!!!!!!!!!!
	print_tensorinfo(&(pair->opt), keeper);
}

void test_blablabla(TensorInfoPair test)
{
	// do nothing
}


// std::ostream& operator << (std::ostream& out, const TensorInfo& tensor) {
// 	printf("--------------------------------------------------------------------------------\n");
//     print_bonds(tensor.get_data());
//     print_couplings(tensor.get_data());
//     printf("\n");
// 	return out;
// }
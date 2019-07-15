// class TwoSiteOverlapCalculator : public OverlapCalculator {
// 	public:
// 		TwoSiteOverlapCalculator
// 	private:
// 		void set_OO_to_contraction(const TensorInfo * A, const TensorInfo * B,
// 				int leg, OverlapObject * OO);
// 		// storage of intermediate results
// 		struct TensorInfo[2] tensmem; // tensor memory
// 		// backup tensor memory for branching tensor calculations
// 		struct TensorInfo[2] tensbmem;
// 		// indices of the last optimization center used
// 		// the indices are stored such that last_optimized
// 		int last_optimized[2];
// };
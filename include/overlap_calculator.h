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

typedef struct T3NSfill {
    struct siteTensor ** data;
    struct bookkeeper * bookie;
} T3NSfill;


#ifdef __cplusplus

// The overlap calculator does not allocate any of the T3NS data operated on.
// It does allocate its own temporary results.
class OverlapCalculator {
    public:
        OverlapCalculator(int test);
        int get_result() { return test; };
    private:
        int test;
};


#else
typedef struct OverlapCalculator OverlapCalculator;
#endif


#ifdef __cplusplus
extern "C" {
#endif

extern void init_overlap_calculator(int test, OverlapCalculator** result);
extern int get_result(OverlapCalculator*);

// #if defined(__STDC__) || defined(__cplusplus)
//   //extern void c_function(Fred*);   /* ANSI C prototypes */
//   extern void init_overlap_calculator(int test, OverlapCalculator** result);
// #else
//   //extern void c_function();        /* K&R style */
//   extern OverlapCalculator** init_overlap_calculator();
// #endif

#ifdef __cplusplus
}
#endif

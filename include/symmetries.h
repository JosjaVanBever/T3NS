/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018-2019 Klaas Gunst <Klaas.Gunst@UGent.be>
    
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
#pragma once

#include "symmetry_z2.h"
#include "symmetry_u1.h"
#include "symmetry_su2.h"
#include "symmetry_pg.h"

/**
 * \file symmetries.h
 * @brief Wrapper for the different symmetries implemented.
 *
 * This is a wrapper for the different symmetries, at this moment we have
 * \f$Z_2, U(1), SU(2),C_1, C_i, C_2, C_s, D_2, C_{2v}, C_{2h}, D_{2h}\f$.
 */

// POINT_GROUP_SYMMETRY is a macro!!!
// If you change this, you should change the symmetrynames in symmetries.c also
enum symmetrygroup { Z2, U1, SU2, POINT_GROUP_SYMMETRY, SENIORITY };

/**
 * @brief Gives the maximal label + 1 of the irreps that can be generated
 * symmetrygroup for a given prop1 and prop2. (prop1 and prop2 should be given
 * so u1 and su2 knows the range)
 *
 * @param [in] prop1 The first array of irreps.
 * @param [in] nr1 The number of irreps in prop1.
 * @param [in] prop2 The second array of irreps.
 * @param [in] nr2 The number of irreps in prop2.
 * @param [in] inc increment between irreps in prop1 and prop2.
 * @param [in] sg The symmetry group of which prop1 and prop2 are irreps.
 * @return returns the maximal label of the irreps that can be generated.
 */
int get_max_irrep(int (*prop1)[MAX_SYMMETRIES], int nr1, 
                  int (*prop2)[MAX_SYMMETRIES], int nr2,
                  enum symmetrygroup sg, int whichsym);

/**
 * @brief Gives the resulting irreps from tensor product of two other irreps
 * belonging to sg.
 *
 * @param [out] min_irrep The lowest label of resulting irrep.
 * @param [out] nr_irreps Number of resulting irreps.
 * @param [out] step Step with which the labels are separated.
 * @param [in] irrep1 The first irrep of the tensorproduct.
 * @param [in] irrep2 The second irrep of the tensorproduct.
 * @param [in] sign -1 if the inverse of irrep2 should be taken, +1 otherwise.
 * @param [in] sg The symmetry group of the irreps.
 */
void tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, int irrep1, 
                    int irrep2, int sign, enum symmetrygroup sg);

/**
 * @brief Returns the string of the symmetrygroup.
 *
 * @param [in] sg The symmetrygroup
 * @return The string of the symmetrygroup.
 */
const char * get_symstring(enum symmetrygroup sg);

void get_allsymstringnames(char * buffer);

/**
 * @brief Searches for an inputted string the right symmetrygroup.
 *
 * @param [in] buffer The string.
 * @param [out] sg The resulting symmetrygroup.
 * @return Returns 1 if successful, otherwise 0.
 */
int which_symmgroup(char * buffer, enum symmetrygroup * sg);


/**
 * @brief Gets the string of the irrep.
 *
 * @param [out] buffer The string.
 * @param [in] sg The symmetrygroup.
 * @param [in] irr The irrep.
 */
void get_irrstring(char * buffer, enum symmetrygroup sg, int irr);

/**
 * @brief Finds the irrep that is inputted in a string.
 *
 * @param [in] buffer The string to read.
 * @param [in] sg The symmetrygroup of which the irrep is an element.
 * @param [out] irr The found irrep.
 * @return Returns 1 if successful, otherwise 0.
 */
int which_irrep(char * buffer, enum symmetrygroup sg, int * irr);

/**
 * @brief Searches for a string in a given array of strings.
 * If the string is just an indexnumber, this is also valid.
 *
 * @param [in] buffer The string to search for.
 * @param [in] arr The array of strings in which to search.
 * @param [in] length The number of elements in the array arr.
 * @param [out] ind The index of which the array is found in the string.
 * @return 1 if the search was successful, 0 otherwise.
 */
int find_str_in_array(char * buffer, const char ** arr, int length, int * ind);

/**
 * @brief The prefactor for when appending a site-operator to a @ref
 * rOperators.
 *
 * @param[in] irrep_arr The different irreps of the symmetries.<br>
 * The irreps are given by
 * \f$[[α', i', β'], [α, i, β] [MPO(α), MPO(i), MPO(β)]]\f$
 *
 * @param[in] is_left The value of rOperators.is_left.
 * @param[in] sgs The symmetrygroups.
 * @param[in] nrsy The number of symmetrygroups.
 * @return The prefactor.
 */
double prefactor_pAppend(const int * (*irrep_arr)[3], int is_left, 
                         const enum symmetrygroup * sgs, int nrsy);

/**
 * @brief The prefactor when making the adjoint of a three-legged T3NS-tensor.
 *
 * For the coupling of the tensor, see @ref siteTensor.
 *
 * @param [in] irreps The different irreps of the symmetries of the current 
 * symmetry sector.<br>
 * The order is given by: \f$α, β, γ\f$.
 * @param [in] c The type of orthogonalized tensor. Options are:
 * * 'c' for an orthogonality center.
 * * '1' for orthogonalized tensors with respect of contraction over β and γ.
 * * '2' for orthogonalized tensors with respect of contraction over α and γ.
 * * '3' for orthogonalized tensors with respect of contraction over α and β.
 * @param [in] sgs The symmetrygroups.
 * @param [in] nrsy The number of symmetrygroups.
 */
double prefactor_adjoint(const int ** irreps, char c, 
                         const enum symmetrygroup * sgs, int nrsy);

/**
 * @brief The prefactor when updating a physical rOperators by contracting with
 * a physical siteTensor and its adjoint.
 *
 * @param [in] irreps The different irreps of the symmetries for the current
 * symmetry sector.<br>
 * The order is the same as the qnumbers order for a physical @ref rOperators.
 * (See table in documentation for @ref rOperators.)
 * @param [in] is_left 1 if updating a left rOperators, otherwise 0.
 * @param [in] sgs The symmetrygroups.
 * @param [in] nrsy The number of symmetrygroups.
 */
double prefactor_pUpdate(const int * (*irrep_arr)[3], int is_left, 
                         const enum symmetrygroup * sgs, int nrsy);

double prefactor_mirror_coupling(int ** irrep_arr, 
                                 const enum symmetrygroup * sgs, int nrsy);

double prefactor_bUpdate(int * (*irrep_arr)[3], int updateCase,
                         const enum symmetrygroup * sgs, int nrsy);

double prefactor_add_P_operator(int * const (*irreps)[3], int isleft, 
                                const enum symmetrygroup * sgs, int nrsy);

double prefactor_combine_MPOs(int * const (*irreps)[3], int * const *irrMPO, 
                              const enum symmetrygroup * sgs, int nrsy,
                              int isdmrg, int extradinge);

/**
 * @brief Calculates the prefactor for a certain permutation of the orbitals of
 * the multi-site tensor.
 *
 * The type of the permutation is given by:
 *      * DMRG:
 *              * 0 : 1 ↔ 2
 *      * T3NS:
 *              * 1 : 2 ↔ 3 (1, 3, 2)
 *              * 2 : 1 ↔ 3 (3, 2, 1)
 *              * 3 : 1 ↔ 2 (2, 1, 3)
 *              * 4 : 1 ↔ 2, 2 ↔ 3 (2, 3, 1)
 *              * 5 : 1 ↔ 2, 1 ↔ 3 (3, 1, 2)
 *
 * The @ref irreps array is given by:
 *      * DMRG:\f$[[α, i, β], [β, j, γ], \mathrm{NULL}, \mathrm{NULL}, [β', i, γ]]\f$<br>
 *      Where the first and the second array are the irreps of the first and
 *      second physical site **after** permutation. The last array are the
 *      irreps of the second physical site **before** permutation.
 *      * T3NS:\f$[[α, i, β], [γ, j, δ], [μ, k, ν], [β, δ, μ], [β', δ', μ']]\f$<br>
 *      Where the first three arrays are the irreps of the first three physical
 *      sites **after** permutation. The fourth array contains the irreps of
 *      the branching site **after** permutation and the last array **before**
 *      permutation.
 *
 * @param [in] irreps The different irreps of bonds.
 * @param [in] permuteType The type of the permutation.
 * @param [in] sgs The symmetrygroups.
 * @param [in] nrsy The number of symmetrygroups.
 * @return The prefactor for this type of permutation.
 */
double prefactor_permutation(int * irreps[5][3], int permuteType, 
                             const enum symmetrygroup * sgs, int nrsy);

/**
 * @brief Returns the prefactor for making the 1-site RDM.
 *
 * In this step, a certain orthocenter is contracted with itself over \f$α\f$
 * and \f$β\f$.
 *
 * The coupling is changed in the following way:
 * \f$(|α〉,|i〉,〈β|),(|β'〉,〈i'|,〈α'|) → (〈i'|, 0, |i〉)\f$
 *
 * **Note :** For graded \f$\mathbb{Z}_2\f$-symmetry, one should need an extra
 * \f$(-1)^{β'}\f$ prefactor. However, this is canceled with the prefactor 
 * needed from the adjoint of a orthonormality center.
 *
 * @param irreps [in] The different irreps of the symmetries of the current 
 * symmetry sector.<br>
 * The order is given by: \f$α, i, β\f$.
 * @param sgs [in] The symmetrygroups.
 * @param nrsy [in] The number of symmetrygroups.
 * @return The prefactor.
 */
double prefactor_1siteRDM(int * (*irreps)[3], const enum symmetrygroup * sgs,
                          int nrsy);


/** @brief The prefactor when initializing an intermediate RDM needed for the
 * calculation of the 2-site RDM's.
 *
 * Two physical tensor are contracted over bond \f$α\f$ or \f$β\f$, the physical 
 * bond is not contracted over.
 *
 * The coupling is changed in the following way:
 * * For contraction over \f$α\f$: \f$(|α〉,|i〉,〈β|),(|β'〉,〈i'|,〈α'|) → 
 *   (|i〉,〈i'|, 〈ii'|), (|ii'〉, 〈β|, |β'〉)\f$
 * * For contraction over \f$β\f$: \f$(|α〉,|i〉,〈β|),(|β'〉,〈i'|,〈α'|) → 
 *   (|i〉,〈i'|, 〈ii'|), (|ii'〉, |α〉, 〈α'|)\f$
 *
 * @param irreps [in] The different irreps of the symmetries of the current 
 * symmetry sector.<br>
 * The order is given by: \f$α, i, β, α', i', β', ii'\f$.
 * @param bond [in] The virtual bond that is left open, thus<br>
 * * 0 if contraction over \f$β\f$
 * * 2 if contraction over \f$α\f$
 * @param sgs [in] The symmetrygroups.
 * @param nrsy [in] The number of symmetrygroups.
 * @return The prefactor.
 */
double prefactor_RDMinterm(int * (*irreps)[7], int bond, 
                           enum symmetrygroup * sgs, int nrsy);

int need_multiplicity(int nrSyms, const enum symmetrygroup * sgs);

int multiplicity(int nrSyms, const enum symmetrygroup * sgs, const int * irreps);

/**
 * @brief returns a buffer of the different symmetry groups inputted.
 *
 * @param [in] sgs The symmetry groups.
 * @param [in] nrSyms The length of the sgs array.
 * @param [in,out] buffer Allocated buffer space.
 */
void get_sgsstring(enum symmetrygroup * sgs, int nrSyms, char * buffer);

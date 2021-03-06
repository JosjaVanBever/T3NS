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

/**
 * \file symmetry_su2.h
 * \brief file for the \f$SU(2)\f$ symmetry.
 *
 * The irreps are labeled as \f$2j\f$. Thus:\n
 * \f$ j = \{0, 1/2, 1, 3/2, 2, \cdots\}\f$ becomes label \f$= \{0, 1, 2, 3, 4, \cdots\}\f$
 */

/**
 * \brief Gives the maximal label + 1 of the irreps that can be generated by SU2.
 *
 * \param [in] prop1 The first array of irreps.
 * \param [in] nr1 The number of irreps in prop1.
 * \param [in] prop2 The second array of irreps.
 * \param [in] nr2 The number of irreps in prop2.
 * \param [in] inc increment between irreps in prop1 and prop2.
 * \return returns the maximal label of the irreps that can be generated.
 */
int SU2_get_max_irrep(int (*prop1)[MAX_SYMMETRIES], int nr1, 
                  int (*prop2)[MAX_SYMMETRIES], int nr2, int whichsym);

/**
 * \brief Gives the resulting irreps from tensor product of two other irreps belonging to sg.
 *
 * \param [out] min_irrep The lowest label of resulting irrep.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [out] step Step with which the labels are separated.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 */
void SU2_tensprod_irrep(int * min_irrep, int * nr_irreps, int * step, 
                        int irrep1, int irrep2);

/**
 * \brief Returns the irrepstring, or INVALID if invalid.
 *
 * \param[out] buffer The string.
 * \param[in] irr The irrep.
 */
void SU2_get_irrstring(char * buffer, int irr);

/**
 * \brief finds which irrep is given in the buffer.
 *
 * \param [in] buffer The buffer.
 * \param [out] irr The irrep.
 * \return 1 if successful, 0 otherwise.
 */
int SU2_which_irrep(char * buffer, int * irr);

double SU2_prefactor_mirror_coupling(const int * symv);

double SU2_prefactor_pAppend(int (*sv)[3], int is_left);

double SU2_prefactor_combine_MPOs(int (*symv)[3], int * symvMPO, int isdmrg, int extradinge);

double SU2_prefactor_bUpdate(int (*symv)[3], int uCase);

/**
 * @brief Returns the prefactor for \f$SU(2)\f$ for making the 1-site RDM.
 *
 * In this step, a certain orthocenter is contracted with itself over \f$α\f$
 * and \f$β\f$.
 *
 * The change of coupling is explained in @ref prefactor_1siteRDM.<br>
 * The prefactor needed when using \f$SU(2)\f$ is \f$\frac{1}{2 j_i + 1}\f$.
 *
 * @param symv [in] The \f$2j\f$ values of the different bonds.<br>
 * The order is given by: \f$α, i, β\f$.
 * @return The prefactor.
 */
double SU2_prefactor_1siteRDM(int * symv);

/** 
 * @brief The prefactor for \f$SU(2)\f$-symmetry when initializing an
 * intermediate RDM needed for the calculation of the 2-site RDM's.
 *
 * The changing of coupling can be found in @ref prefactor_RDMinterm.<br>
 * The prefactors arising from this change of coupling are:
 * * For contraction over \f$α: (-1)^{α + í' - β'} [ii']^2 [β][β']
 *   \begin{Bmatrix}
 *   i & i' & ii'\\
 *   β' & β & α
 *   \end{Bmatrix}\f$
 * * For contraction over \f$β: (-1)^{ii'  - β' - α' - i} [ii']^2 [α][α']
 *   \begin{Bmatrix}
 *   i & i' & ii'\\
 *   α' & α & β
 *   \end{Bmatrix}\f$
 *
 * @param symvalues [in] The parities of the different bonds.<br>
 * The order is given by: \f$α, i, β, α', i', β', ii'\f$.
 * @param bond [in] The virtual bond that is left open, thus<br>
 * * 0 if contraction over \f$β\f$
 * * 2 if contraction over \f$α\f$
 * @return The prefactor.
 */
double SU2_prefactor_RDMinterm(int * symvalues, int bond);

int SU2_multiplicity(const int irrep);

/**
 * @brief Calculates the prefactor for a certain permutation of the orbitals of
 * the multi-site tensor.
 *
 * The type of the permutation and their prefactor is given by:
 *      * DMRG:
 *              * 0 : 1 ↔ 2:<br>
 *              \f$(-1)^{j_i + j_j - (j_β + j_{β'})}[j_β][j_{β'}]
 *              \begin{Bmatrix}
 *              j_i & j_γ & j_{β'}\\
 *              j_j & j_α & j_β
 *              \end{Bmatrix}
 *              \f$
 *      * T3NS:
 *              * 1 : 2 ↔ 3 (1, 3, 2):<br>
 *              \f$(-1)^{2(j_δ - j_{δ'})}[j_δ][j_{δ'}][j_μ][j_{μ'}]
 *              \begin{Bmatrix}
 *              j_γ & j_k & j_{δ'}\\
 *              j_j & j_ν & j_{μ'}\\
 *              j_δ & j_μ & j_β
 *              \end{Bmatrix}
 *              \f$
 *
 *              * 2 : 1 ↔ 3 (3, 2, 1):<br>
 *              \f$(-1)^{(j_β - j_{β'}) - (j_μ - j_{μ'})}[j_β][j_{β'}][j_μ][j_{μ'}]
 *              \begin{Bmatrix}
 *              j_α & j_k & j_{β'}\\
 *              j_i & j_ν & j_{μ'}\\
 *              j_β & j_μ & j_δ
 *              \end{Bmatrix}
 *              \f$
 *
 *              * 3 : 1 ↔ 2 (2, 1, 3):<br>
 *              \f$(-1)^{(j_j - j_i) + (j_δ - j_{δ'})}[j_β][j_{β'}][j_δ][j_{δ'}]
 *              \begin{Bmatrix}
 *              j_α & j_j & j_{β'}\\
 *              j_i & j_γ  & j_{δ'}\\
 *              j_β & j_δ & j_μ
 *              \end{Bmatrix}
 *              \f$
 *
 *              * 4 : 1 ↔ 2, 2 ↔ 3 (2, 3, 1):<br>
 *              \f$ Σ_J
 *              (-1)^{(j_k - j_i) + (J - j_{δ'})}[j_β][j_{β'}][J][j_{δ'}]
 *              \begin{Bmatrix}
 *              j_α & j_k & j_{β'}\\
 *              j_i & j_γ  & j_{δ'}\\
 *              j_β & J & j_{μ'}
 *              \end{Bmatrix}
 *              (-1)^{2(j_δ - J)}[j_δ][J][j_μ][j_{μ'}]
 *              \begin{Bmatrix}
 *              j_γ & j_k & J\\
 *              j_j & j_ν & j_{μ'}\\
 *              j_δ & j_μ & j_β
 *              \end{Bmatrix}
 *              \f$
 *
 *              * 5 : 1 ↔ 2, 1 ↔ 3 (3, 1, 2):<br>
 *              \f$ Σ_J
 *              (-1)^{(j_j - j_k) + (j_δ - j_{δ'})}[J][j_{β'}][j_δ][j_{δ'}]
 *              \begin{Bmatrix}
 *              j_α & j_j & j_{β'}\\
 *              j_k & j_γ  & j_{δ'}\\
 *              J & j_δ & j_{μ'}
 *              \end{Bmatrix}
 *              (-1)^{(j_β - J) - (j_μ - j_{μ'})}[j_β][J][j_μ][j_{μ'}]
 *              \begin{Bmatrix}
 *              j_α & j_k & J\\
 *              j_i & j_ν & j_{μ'}\\
 *              j_β & j_μ & j_δ
 *              \end{Bmatrix}
 *              \f$
 *
 * For the order of the symv, see the documentation on prefactor_permutation().
 *
 * @param [in] symv The different \f$2j\f$ values of the bonds.
 * @param [in] permuteType The type of the permutation.
 * @return The prefactor for this type of permutation.
 */
double SU2_prefactor_permutation(int (*symv)[3], int permuteType);

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
 * \file symmetry_pg.h
 * \brief Header file for the abelian point groups implemented.
 *
 *  This header contains the symmetry group and irrep conventions.
 *  The program requires Abelian point groups with real character tables,
 *  with hence \f$I_{\alpha} \otimes I_{\alpha} = I_{trivial}\f$.
 *
 *  \section irreps_conv Irrep conventions
 *  
 *  The same conventions as in Psi4 (beta5) are used. For convenience, they are listed below:\n
 *  \latexonly
 *  \begin{tabular}{|l|cccccccc|}
 *  \hline
 *    Symmetry Conventions & \multicolumn{8}{c|}{ Irrep Number \& Name } \\
 *  \hline
 *    Group Number \& Name & 0  & 1   & 2   & 3   & 4  & 5   & 6   & 7   \\
 *  \hline
 *    0: c1                & A  &     &     &     &    &     &     &     \\
 *    1: ci                & Ag & Au  &     &     &    &     &     &     \\
 *    2: c2                & A  & B   &     &     &    &     &     &     \\
 *    3: cs                & A' & A'' &     &     &    &     &     &     \\
 *    4: d2                & A  & B1  & B2  & B3  &    &     &     &     \\
 *    5: c2v               & A1 & A2  & B1  & B2  &    &     &     &     \\
 *    6: c2h               & Ag & Bg  & Au  & Bu  &    &     &     &     \\
 *    7: d2h               & Ag & B1g & B2g & B3g & Au & B1u & B2u & B3u \\
 *  \hline
 *  \end{tabular}
 *  \endlatexonly
 *  \htmlonly
 *  <table border="1">
 *  <tr><td> Symmetry Conventions </td><td colspan="8"> Irrep Number & Name </td></tr>
 *  <tr><td> Group Number & Name </td><td> 0  </td><td> 1   </td><td> 2   </td><td> 3   </td><td> 4  </td><td> 5   </td><td> 6   </td><td> 7   </td></tr>
 *  <tr><td> 0: c1               </td><td> A  </td><td>     </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
 *  <tr><td> 1: ci               </td><td> Ag </td><td> Au  </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
 *  <tr><td> 2: c2               </td><td> A  </td><td> B   </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
 *  <tr><td> 3: cs               </td><td> A' </td><td> A'' </td><td>     </td><td>     </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
 *  <tr><td> 4: d2               </td><td> A  </td><td> B1  </td><td> B2  </td><td> B3  </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
 *  <tr><td> 5: c2v              </td><td> A1 </td><td> A2  </td><td> B1  </td><td> B2  </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
 *  <tr><td> 6: c2h              </td><td> Ag </td><td> Bg  </td><td> Au  </td><td> Bu  </td><td>    </td><td>     </td><td>     </td><td>     </td></tr>
 *  <tr><td> 7: d2h              </td><td> Ag </td><td> B1g </td><td> B2g </td><td> B3g </td><td> Au </td><td> B1u </td><td> B2u </td><td> B3u </td></tr>
 *  </table>
 *  \endhtmlonly
 *  Note that these conventions allow to use the XOR operation for irrep multiplication.
 */
/* NOTE : maybe I just can leave out C1 since this is a trivial symmetry */
#define POINT_GROUP_SYMMETRY C1, Ci, C2, Cs, D2, C2v, C2h, D2h

/**
 * \brief Gives the maximal label + 1 of the irreps that can be generated by the point group.
 *
 * \param [in] pg The point group symmetry.
 * \return returns the maximal label of the irreps that can be generated.
 */
int PG_get_max_irrep(int pg);

/**
 * \brief Gives the resulting irreps from tensor product of two other irreps belonging to sg.
 *
 * \param [out] min_irrep The lowest label of resulting irrep.
 * \param [out] nr_irreps Number of resulting irreps.
 * \param [out] step Step with which the labels are separated.
 * \param [in] irrep1 The first irrep of the tensorproduct.
 * \param [in] irrep2 The second irrep of the tensorproduct.
 */
void PG_tensprod_irrep(int *min_irrep, int *nr_irreps, int *step, 
                       int irrep1, int irrep2);

/**
 * \brief Returns the irrepstring, or INVALID if invalid.
 *
 * \param [out] buffer The resulting string. 
 * \param [in] pg The point group symmetry.
 * \param [in] irr The irrep.
 */
void PG_get_irrstring(char * buffer, int pg, int irr);

/**
 * \brief finds which irrep is given in the buffer.
 *
 * \param [in] buffer The buffer.
 * \param [in] pg The point group symmetry.
 * \param [out] irr The irrep.
 * \return 1 if successful, 0 otherwise.
 */
int PG_which_irrep(char * buffer, int pg, int *irr);

int fcidump_to_psi4(const int fcidumpirrep, const int pg_symm);

/*
    T3NS: an implementation of the Three-Legged Tree Tensor Network algorithm
    Copyright (C) 2018 Klaas Gunst <Klaas.Gunst@UGent.be>
    
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
#include <stdlib.h>
#include <stdio.h>

#include "rOperators.h"
#include "debug.h"
#include "network.h"
#include "instructions.h"
#include "hamiltonian.h"
#include "lapack.h"
#include "sort.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void unique_rOperators_append_phys(struct rOperators * const uniquerops,
     const int * const instructions, const int * const hamsymsecs_of_new, const int nr_instructions,
     const struct rOperators * const fullrops);

static void initialize_unique_rOperators_append_phys(struct rOperators * const uniquerops, 
    const int bond_of_operator, const int is_left, const int * const instructions, 
    const int * const hamsymsecs_of_new, const int nr_instructions);

static double calculate_prefactor_append_physical(const int indexes[], const struct symsecs 
    symarr[], const int is_left);

static void get_separate_qnumbers(int indexes[], const QN_TYPE * const qnumbers, int maxdims[], 
    const int couplings);

static int get_compressed_instructions(const int * const instructions, const int nr_instructions, 
    const int * const prev_hss, const int site, const int * const new_hss, int ** const compr_instr,
    int ** const compr_hss);

static QN_TYPE change_new_to_old_ops(const QN_TYPE new, const int dimold, const int dimnew, 
    int * const oldtonew);

static QN_TYPE change_old_to_new_site(const QN_TYPE old, const int dimold, const int maxdims[],
    int * const oldtonew, const int is_left);

static int * make_oldtonew(const struct symsecs * const internalss, const int bond);

#ifdef DEBUG
static int consistency_check_for_update_physical(const struct rOperators * const rops, 
    const struct siteTensor * const tens);
#endif

/* ========================================================================== */

void rOperators_append_phys(struct rOperators * const newrops, const struct rOperators * 
    const oldrops)
{
  /**
   * instructions are given in a nr_instructions * 3 array : ops_index, site_op_index, result_index.
   * in the prefactor array, the prefactor of every instruction is stored, in the hamsymsecs_of_new
   * array the index of the irreps of the new renormalized operators are given in.
   *
   * This is needed because for only U(1) there was no ambiguity to which irrep the old_op + site_op
   * should couple to, however this is not the case for SU(2).
   *
   * In the oldrops: some hermitian conjugates are not yet included which we need for the making
   * of newrops, so first we need to expand them.
   *
   * The instructions make a condensed set of renormalized operators.
   */
  int *instructions, *hamsymsecs_of_new, nr_instructions;
  double *prefactors;

  /* Here the result of all the unique physical appends is stored */
  struct rOperators uniquerops; 

  /* hamsymsecs_of_new is an extra array with the hamsymsecs index of every resulting operator */
  fetch_pUpdate(&instructions, &prefactors, &hamsymsecs_of_new, &nr_instructions, 
      oldrops->bond_of_operator, oldrops->is_left);
  unique_rOperators_append_phys(&uniquerops, instructions, hamsymsecs_of_new, 
      nr_instructions, oldrops);

  sum_unique_rOperators(newrops, &uniquerops, instructions, hamsymsecs_of_new, prefactors, 
      nr_instructions);
  destroy_rOperators(&uniquerops);

  safe_free(instructions);
  safe_free(prefactors);
  safe_free(hamsymsecs_of_new);
}

void update_rOperators_physical(struct rOperators * const rops, const struct siteTensor * 
    const tens, const struct symsecs * const internalss)
{
  /* THE SYMMSEC OF THE INTERNAL BOND OF ROPS SHOULD ALREADY BEEN UPDATED TO AN EXTERNAL SYMMSEC.
   * This means not all dims are necessarily equal to 1 (some may be even zero).
   *
   * So for this function, we transform the renormalized operators from a 7-indexed operator
   * to a 3-indexed operator.
   *
   * A certain prefactor should be included for this transform.
   *
   * The hermitian of the tensor is not explicitely constructed, however, the prefactor linked
   * to making it hermitian is.
   *
   * We have the following couplings for the renormalized ops, for tens and for tens_hermitian:
   *
   * renormalized ops (LEFT): (A asterisk means it is an out-bond)
   *   begin:
   *      ---[bra(alpha), bra(i) , bra(beta)* ,
   *           bra(beta) , MPO*   , ket(beta)* ,
   *           ket(beta) , ket(i)*, ket(alpha)*]
   *   end:
   *      ---[bra*(beta), MPO, ket(beta)]
   *
   * renormalized ops (RIGHT):
   *   begin:
   *      ---[bra(alpha) , bra(i) , bra(beta)*,
   *           bra(alpha)*, MPO*   , ket(alpha),
   *           ket(beta)  , ket(i)*, ket(alpha)*]
   *   end:
   *      ---[bra(alpha), MPO, ket*(alpha)]
   *
   * tens:
   *      ---[ket(alpha), ket(i), ket(beta)*]
   * tens_hermitian:
   *      ---[bra(beta), bra(i)*, bra(alpha)*]
   */
  const int oldcoupling = 3;
  const int is_left = rops->is_left;
  int ** tmpblockbegin;
  struct rOperators updated_rops;
  struct symsecs symarr[4];
  int i, new_sb;
  int old_sb = 0;

  /* Do for left renormalized operators  : tens_herm_sb.T x old_sb x tens_sb
   * Do for right renormalized operators : tens_herm_sb x old_sb x tens_sb.T */
  const char TRANS1 = is_left ? 'T' : 'N';
  const char TRANS2 = is_left ? 'N' : 'T';
  const char NOTRANS = 'N';
  const double ONE = 1;
  const double ZERO = 0;

  const int dimoldint = internalss->nrSecs; 
  int * oldtonew = make_oldtonew(internalss, rops->bond_of_operator);
  int maxdims[3];
  int indices[3];
  siteTensor_give_indices(tens, indices);
  get_symsecs_arr(3, symarr, indices);
  get_symsecs(&symarr[3], get_hamiltonianbond(rops->bond_of_operator));

  /* Check if the bonds of the renormalized operators and the bonds of the tensor are 
   * corresponding, put here an assert that checks all bonds of tens and rops if they have the
   * right correspondence and all other weird as shit. ten should be a 1 site tensor and so on ...*/
  assert(consistency_check_for_update_physical(rops, tens));

  for (i = 0; i < 3; ++i) maxdims[i] = symarr[i].nrSecs;

  /* initialize the three-indexed renormalized operator */
  init_rOperators(&updated_rops, &tmpblockbegin, rops->bond_of_operator, is_left, 0);
  updated_rops.nrops      = rops->nrops;
  updated_rops.hss_of_ops = safe_malloc(updated_rops.nrops, int);
  updated_rops.operators  = safe_malloc(updated_rops.nrops, struct sparseblocks);
  for (i = 0; i < rops->nrops; ++i)
  {
    const int currhss = rops->hss_of_ops[i];
    const int nr_blocks = rOperators_give_nr_blocks_for_hss(&updated_rops, currhss);
    struct sparseblocks * const currBlock = &updated_rops.operators[i];
    int j;
    updated_rops.hss_of_ops[i] = currhss;

    currBlock->beginblock = safe_malloc(nr_blocks + 1, int);

    for (j = 0; j < nr_blocks + 1; ++j)
      currBlock->beginblock[j] = tmpblockbegin[currhss][j];
    currBlock->tel = safe_calloc(currBlock->beginblock[nr_blocks], EL_TYPE);
  }
  for (i = 0; i < updated_rops.nrhss; ++i)
    safe_free(tmpblockbegin[i]);
  safe_free(tmpblockbegin);

  /* Loop over the different symmetryblocks of the new renormalized operators. */
  for (new_sb = 0; new_sb < updated_rops.begin_blocks_of_hss[updated_rops.nrhss]; ++new_sb)
  {
    int newhss = 0;
    const QN_TYPE newqnumber = updated_rops.qnumbers[new_sb];
    const QN_TYPE newtooldqn = change_new_to_old_ops(newqnumber, dimoldint, maxdims[2 * is_left],
        oldtonew);

    /* Search the hamiltonian_symsec of the current symmetryblock */
    while (new_sb >= updated_rops.begin_blocks_of_hss[newhss + 1]) ++newhss;

    /* NOTE: This can be probably be done more efficient, by not looping over all old_sb since
     * the ones that pass the next if statement are all consecutive */
    for (; old_sb < rops->begin_blocks_of_hss[rops->nrhss]; ++old_sb)
    {
      /* if newqnumber corresponds with the last qnumber of old_sb than you should do something */
      //if (newqnumber == rops->qnumbers[oldcoupling * old_sb + oldcoupling - 1])
      if (newtooldqn == rops->qnumbers[oldcoupling * old_sb + oldcoupling - 1])
      {
        /* T3NS qnumber is the second qnumber of old_sb, T3NS herm qnumber is the first one */
        const QN_TYPE T3NSoldqnumber     = rops->qnumbers[oldcoupling * old_sb + 1];
        const QN_TYPE T3NSoldhermqnumber = rops->qnumbers[oldcoupling * old_sb];
        const QN_TYPE T3NSqnumber     = change_old_to_new_site(T3NSoldqnumber, dimoldint, 
            maxdims, oldtonew, is_left);
        const QN_TYPE T3NShermqnumber = change_old_to_new_site(T3NSoldhermqnumber, dimoldint, 
            maxdims, oldtonew, is_left);
        double prefactor;
        const int tens_sb       = siteTensor_search_qnumber(T3NSqnumber, tens);
        const int tens_herm_sb  = siteTensor_search_qnumber(T3NShermqnumber, tens);
        const int new_sb_tens   = new_sb - updated_rops.begin_blocks_of_hss[newhss];
        const int old_sb_tens   = old_sb - rops->begin_blocks_of_hss[newhss];

        const struct sparseblocks * prev_operator = &rops->operators[0];
        struct sparseblocks * curr_operator       = &updated_rops.operators[0];

        /* array consisting of indexes of:
         *    bra(alpha) bra(i) bra(beta) ket(alpha) ket(i) ket(beta) MPO */
        const int ind_arr[7] = 
        { T3NShermqnumber % maxdims[0], 
          (T3NShermqnumber / maxdims[0]) % maxdims[1],
          (T3NShermqnumber / maxdims[0]) / maxdims[1],
          T3NSqnumber % maxdims[0], 
          (T3NSqnumber / maxdims[0]) % maxdims[1],
          (T3NSqnumber / maxdims[0]) / maxdims[1],
          newhss };
        const int * irrep_arr[7] = 
        { symarr[0].irreps[ind_arr[0]],
          symarr[1].irreps[ind_arr[1]],
          symarr[2].irreps[ind_arr[2]],
          symarr[0].irreps[ind_arr[3]],
          symarr[1].irreps[ind_arr[4]],
          symarr[2].irreps[ind_arr[5]],
          symarr[3].irreps[ind_arr[6]] };

        /* new_sb is size of NxN'
         * old_sb is size of MxM'
         * tens_sb is size of M'xN' and tens_herm_sb is MxN for left renormalized operators
         * tens_sb is size of N'xM' and tens_herm_sb is NxM for right renormalized operators */
        const int N  = symarr[2 * is_left].dims[ind_arr[2 * is_left]];
        const int N2 = symarr[2 * is_left].dims[ind_arr[3 + 2 * is_left]];
        const int M  = symarr[2 * !is_left].dims[ind_arr[2 * !is_left]] *
          symarr[1].dims[ind_arr[1]];
        const int M2 = symarr[2 * !is_left].dims[ind_arr[3 + 2 * !is_left]] *
          symarr[1].dims[ind_arr[4]];

        /* Now calculate the most efficient way to execute the two dgemms.
         * For both left and right rops:  NxMxM' + NxM'xN' or MxM'xN' + NxMxN' */
        const int dgemm_order = N * M2 * (M + N2) > M * N2 * (M2 + N);
        double * const tens_block      = get_tel_block(&tens->blocks, tens_sb);
        double * const tens_herm_block = get_tel_block(&tens->blocks, tens_herm_sb);
        double * workmem;

        /* skip to next symmetryblock if tens_sb or tens_herm_sb is empty */
        if (get_size_block(&tens->blocks, tens_sb) == 0 || 
            get_size_block(&tens->blocks, tens_herm_sb) == 0)
          continue;

        assert(get_size_block(&tens->blocks, tens_sb) == N2 * M2 ||
            get_size_block(&tens->blocks, tens_herm_sb) == N * M);

        workmem = safe_malloc(dgemm_order ? M * N2 : N * M2, double);
        prefactor  = prefactor_adjoint(irrep_arr, is_left ? '3' : '1', bookie.sgs, 
            bookie.nrSyms);
        prefactor *= prefactor_pUpdate(irrep_arr, is_left, bookie.sgs, 
            bookie.nrSyms);

        /* Now the intensive part happens...
         * 
         * Do for left renormalized operators  : tens_herm_sb.T x old_sb x tens_sb
         * Do for right renormalized operators : tens_herm_sb x old_sb x tens_sb.T */
        for (i = 0; i < rops->nrops; ++i, ++prev_operator, ++curr_operator)
        {
          double * old_block;
          double * new_block;
          if (rops->hss_of_ops[i] != newhss)
            continue;

          old_block = get_tel_block(prev_operator, old_sb_tens);
          new_block = get_tel_block(curr_operator, new_sb_tens);
          /* the symblock of the old and new operator can be skipped if empty */
          if (get_size_block(prev_operator, old_sb_tens) == 0 ||
             get_size_block(curr_operator, new_sb_tens) == 0)
          continue;

          /* the symblock of the old and new operator should match with the size from the symsecs */
          assert(get_size_block(prev_operator, old_sb_tens) == M * M2 ||
              get_size_block(curr_operator, new_sb_tens) == N * N2);

          /* First way of contracting : old_sb x tens_sb.TRANS2 --> tens_herm_sb.TRANS1 x workmem */
          if (dgemm_order)
          {
            dgemm_(&NOTRANS, &TRANS2, &M, &N2, &M2, &ONE, old_block, &M, tens_block, 
                TRANS2 == 'N' ? &M2 : &N2, &ZERO, workmem, &M);
            dgemm_(&TRANS1, &NOTRANS, &N, &N2, &M, &prefactor, tens_herm_block, TRANS1 == 'N' ? &N:
                &M, workmem, &M, &ONE, new_block, &N);
          }
          /* Second way of contracting : tens_herm_sb.TRANS1 x old_sb --> workmem x tens_sb.TRANS2*/
          else
          {
            dgemm_(&TRANS1, &NOTRANS, &N, &M2, &M, &ONE, tens_herm_block, TRANS1 == 'N' ? &N : &M,
                old_block, &M, &ZERO, workmem, &N);
            dgemm_(&NOTRANS, &TRANS2, &N, &N2, &M2, &prefactor, workmem, &N, tens_block, 
                TRANS2 == 'N' ? &M2 : &N2, &ONE, new_block, &N);
          }
        }
        safe_free(workmem);
      }
      else if (newtooldqn < rops->qnumbers[oldcoupling * old_sb + oldcoupling - 1])
        break;
    }
  }

  clean_symsecs_arr(3, symarr, indices);
  clean_symsecs(&symarr[3], get_hamiltonianbond(rops->bond_of_operator));
  safe_free(oldtonew);
  destroy_rOperators(rops);
  *rops = updated_rops;
  for (i = 0; i < rops->nrops; ++i)
    kick_zero_blocks(&rops->operators[i], rOperators_give_nr_blocks_for_operator(rops, i));
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void unique_rOperators_append_phys(struct rOperators * const uniquerops,
    const int * const instructions, const int * const hamsymsecs_of_new, const int nr_instructions,
     const struct rOperators * const prevrops)
{
  /* The Lsite for right, and Rsite for left */
  const int site = netw.bonds[prevrops->bond_of_operator][prevrops->is_left];
  int bonds[3];
  int tmpbond[9];
  struct symsecs symarr[9];
  int newsec;
  const int couplings = 3;
  int maxdims[3 * couplings];
  int qnbonds[3 * couplings];
  int *compr_instr;
  int *compr_hss;
  int nrinstr;
  get_bonds_of_site(site,  bonds);

  /* bra(alpha), bra(i), bra(beta) */
  tmpbond[0] = get_braT3NSbond(bonds[0]);
  tmpbond[1] = get_braT3NSbond(bonds[1]);
  tmpbond[2] = get_braT3NSbond(bonds[2]);

  /* ket(alpha), ket(i), ket(beta) */
  tmpbond[3] = get_ketT3NSbond(bonds[0]);
  tmpbond[4] = get_ketT3NSbond(bonds[1]);
  tmpbond[5] = get_ketT3NSbond(bonds[2]);

  /* hamsymsec_alpha, hamsymsec_site, hamsymsec_beta */
  tmpbond[6] = get_hamiltonianbond(bonds[0]);
  tmpbond[7] = get_hamiltonianbond(bonds[1]);
  tmpbond[8] = get_hamiltonianbond(bonds[2]);

  get_symsecs_arr(9, symarr, tmpbond);
 
  assert(is_psite(site));

  initialize_unique_rOperators_append_phys(uniquerops, bonds[2 * prevrops->is_left], 
      prevrops->is_left, instructions, hamsymsecs_of_new, nr_instructions);

  assert(rOperators_give_nr_of_couplings(prevrops) == 1 && 
      rOperators_give_nr_of_couplings(uniquerops) == couplings);

  rOperators_give_qnumberbonds(uniquerops, qnbonds);
  get_maxdims_of_bonds(maxdims, qnbonds, 3 * couplings);
  nrinstr = get_compressed_instructions(instructions, nr_instructions, prevrops->hss_of_ops, site, 
      hamsymsecs_of_new, &compr_instr, &compr_hss);

  /* Now loop over different symsecs of uniquerops, find the symsecs of the original that correspond
   * with it and of the siteoperator, and loop over all these possibilities also.
   * After this, calculate prefactor for transform for these sectors.
   * Now loop over the different instructions, and see which need these trasforms...
   *
   * First loop over all the new symsecs, newsec is label number according to rops */
  for (newsec = 0; newsec < uniquerops->begin_blocks_of_hss[uniquerops->nrhss]; ++newsec)
  {
    /* The indexes of the symmsecs of all 9 bonds involved!
     * Column major stored.
     * order is [bra(alpha), ket(alpha), hamsymsec_alpha,
     *            bra(i)  , ket(i)  , hamsymsec_site ,
     *            bra(beta) , ket(beta) , hamsymsec_beta]
     * Where Hamsymsec_alpha is hamsymsec_old for left operator and hamsymsec_new for right operator
     * Where Hamsymsec_beta  is hamsymsec_new for left operator and hamsymsec_old for right operator
     */
    int indexes[9];
    int nr_of_prods;
    int prod;
    int *possible_prods;
    int hamsymsec_new;

    QN_TYPE * const newqnumber = &uniquerops->qnumbers[newsec * couplings];

    /* This finds the qnumber passed in indices of the different bonds.
     * Due to the order of bonds in the renormalized ops, the indices are immediately stored
     * in the right order in the indexes array for the first two couplings. */
    get_separate_qnumbers(indexes, newqnumber, maxdims, couplings);

    /* Hamsymsec_new is Hamsymsec_beta for left, Hamsymsec_new is Hamsymsec_alpha for right.
     * And is stored in the last index element by the get_separate_qnumbers function call. */
    hamsymsec_new = indexes[8];
    indexes[6 + 2 * prevrops->is_left] = hamsymsec_new;

    /* This function should decide which Hamsymsec_site and Hamsymsec_old I need to make the 
     * Hamsymsecnew. This can be pretty much hardcoded sinds this is only hameltonianspecific.
     * The site is also passed to say that the second bond is a site operator bond, cuz I will need 
     * the same function for operatormerging at a branching tensor */
    tprods_ham(&nr_of_prods, &possible_prods, hamsymsec_new, site);

    for (prod = 0; prod < nr_of_prods; ++prod)
    {
      double prefactor;
      const int hamsymsec_site = possible_prods[prod * 2];
      const int hamsymsec_old  = possible_prods[prod * 2 + 1];
      int instr;
      int * curr_instr = compr_instr;
      int * curr_hss = compr_hss;

      /* oldqnumber is bra(alpha), ket(alpha), MPO for left
       * and           bra(beta), ket(beta), MPO for right */ 
      const QN_TYPE oldqnumber = indexes[2 * !prevrops->is_left] + 
        indexes[3 + 2 * !prevrops->is_left] * maxdims[2 * !prevrops->is_left] +
        hamsymsec_old * maxdims[2 * !prevrops->is_left] * maxdims[3 + 2 * !prevrops->is_left];
      assert(maxdims[2 * !prevrops->is_left] == maxdims[3 + 2 * !prevrops->is_left]);

      /* the index of the old and the new block */
      int oldblock;
      int newblock;

      /* pointer to the first operator in the row */
      struct sparseblocks * nextBlock = &uniquerops->operators[0]; 
      /* indexes is completely filled in now */
      indexes[6 + 2 * !prevrops->is_left] = hamsymsec_old;
      indexes[7]                          = hamsymsec_site;

      /* This function calculates the prefactor for this symsec manipulation, for all symmetries */
      prefactor = calculate_prefactor_append_physical(indexes, symarr, uniquerops->is_left);

      oldblock = qnbsearch(&oldqnumber, 1,
          rOperators_give_qnumbers_for_hss(prevrops, hamsymsec_old), 
          rOperators_give_nr_of_couplings(prevrops),
          rOperators_give_nr_blocks_for_hss(prevrops, hamsymsec_old));

      /* symsec not found */
      if (oldblock == -1 || COMPARE_ELEMENT_TO_ZERO(prefactor))
        continue;
      newblock = newsec - uniquerops->begin_blocks_of_hss[hamsymsec_new];

      /* Now loop over the different instructions and only execute the unique ones. */
      for (instr = 0; instr < nrinstr; ++instr, ++nextBlock, curr_instr += 3, curr_hss += 3)
      {
        const int prevoperator = curr_instr[0];
        const int siteoperator = curr_instr[1];
#ifdef DEBUG
        const int nextoperator = curr_instr[2];
#endif
        const struct sparseblocks * const prevBlock = &prevrops->operators[prevoperator];

        assert(prevoperator < prevrops->nrops);

        assert(hamsymsecs_of_new[nextoperator] == 
            uniquerops->hss_of_ops[nextBlock - uniquerops->operators]);
        /* If hamsymsec_old does not correspond with the hamsymsec of the previous operator passed
         * or hamsymsec_site does not correspond with the hamsymsec of the site operator passed
         * or hamsymsec_new does not correspond with the hamsymsec of the new operator passed
         * then you can just skip this instruction, because the relevant symsec manipulation does
         * not occur in this one!
         *
         * If not we can start appending */
        if (curr_hss[0] == hamsymsec_old && curr_hss[1] == hamsymsec_site &&
            curr_hss[2] == hamsymsec_new)
        {
          const int N = get_size_block(prevBlock, oldblock);

          /* This function gets the bra(i), ket(i) element of siteoperator
           * in symsec specified by indexes[1] and indexes[4] */
          const double site_el = prefactor * el_siteop(siteoperator, indexes[1], 
              indexes[4]);
          int j;
          EL_TYPE * const prevTel = get_tel_block(prevBlock, oldblock);
          EL_TYPE * const nextTel = get_tel_block(nextBlock, newblock);

          assert(N == 0 || N == get_size_block(nextBlock, newblock));

          for (j = 0; j < N; ++j) nextTel[j] = site_el * prevTel[j];
        }
      }

      /* check if i looped over all the uniqueoperators */
      assert(nextBlock - uniquerops->operators == uniquerops->nrops);
    }
    safe_free(possible_prods);
  }

  clean_symsecs_arr(9, symarr, tmpbond);
  safe_free(compr_instr);
  safe_free(compr_hss);
}

static void initialize_unique_rOperators_append_phys(struct rOperators * const uniquerops, 
    const int bond_of_operator, const int is_left, const int * const instructions, 
    const int * const hamsymsecs_of_new, const int nr_instructions)
{
  int i;
  int count;
  int **nkappa_begin_temp;

  init_rOperators(uniquerops, &nkappa_begin_temp, bond_of_operator, is_left, 1);

  /* counting number of uniquerops */
  count = 0;
  for (i = 0; i < nr_instructions; ++i)
  {
    const int prevoperator = instructions[i * 3 + 0];
    const int siteoperator = instructions[i * 3 + 1];
    const int nextoperator = instructions[i * 3 + 2];
    if (i == 0 || prevoperator != instructions[(i - 1) * 3 + 0] ||
        siteoperator != instructions[(i - 1) * 3 + 1] || 
        hamsymsecs_of_new[nextoperator] != hamsymsecs_of_new[instructions[(i-1) * 3 + 2]])
      ++count;
  }
  uniquerops->nrops = count;

  /* initializing the hamsymsecs */
  uniquerops->hss_of_ops = safe_malloc(uniquerops->nrops, int);
  count = 0;
  for (i = 0; i < nr_instructions; ++i)
  {
    const int prevoperator = instructions[i * 3 + 0];
    const int siteoperator = instructions[i * 3 + 1];
    const int nextoperator = instructions[i * 3 + 2];
    if (i == 0 || prevoperator != instructions[(i - 1) * 3 + 0] ||
        siteoperator != instructions[(i - 1) * 3 + 1] || 
        hamsymsecs_of_new[nextoperator] != hamsymsecs_of_new[instructions[(i-1) * 3 + 2]])
      uniquerops->hss_of_ops[count++] = hamsymsecs_of_new[nextoperator];
  }
  assert(count == uniquerops->nrops);

  /* initializing the stensors */
  uniquerops->operators = safe_malloc(uniquerops->nrops, struct sparseblocks);
  for (i = 0; i < uniquerops->nrops; ++i)
  {
    /* The current operator and current hamsymsec */
    struct sparseblocks * const blocks = &uniquerops->operators[i];
    const int currhss                  = uniquerops->hss_of_ops[i]; 
    const int N                        = rOperators_give_nr_blocks_for_hss(uniquerops, currhss);

    init_sparseblocks(blocks, nkappa_begin_temp[currhss], N, 'c');
  }

  for (i = 0; i < uniquerops->nrhss; ++i)
    safe_free(nkappa_begin_temp[i]);
  safe_free(nkappa_begin_temp);
}

static double calculate_prefactor_append_physical(const int indexes[], const struct symsecs 
    symarr[], const int is_left)
{
  /* The indexes of the symmsecs of all 9 bonds involved!
   * Column major stored.
   * order is [bra(alpha), ket(alpha), hamsymsec_alpha,
   *            bra(i)  , ket(i)  , hamsymsec_site ,
   *            bra(beta) , ket(beta) , hamsymsec_beta]
   * Where Hamsymsec_alpha is hamsymsec_old for left operator and hamsymsec_new for right operator
   * Where Hamsymsec_beta  is hamsymsec_new for left operator and hamsymsec_old for right operator
   *
   * As can be seen in my notes for example at page 13, for the wigner9j */
  double prefactor = 1;
  int symvalues[9];
  int i, j;

 for (i = 0; i < bookie.nrSyms; ++i)
  {
    for (j = 0; j < 9; ++j)
    {
      symvalues[j] = symarr[j].irreps[indexes[j]][i];
    }
    prefactor *= prefactor_pAppend(symvalues, is_left, bookie.sgs[i]);
  }

  return prefactor;
}

static void get_separate_qnumbers(int indexes[], const QN_TYPE * const qnumbers, int maxdims[], 
    const int couplings)
{
  int i, j;
  int count = 0;

  for (i = 0; i < couplings; ++i)
  {
    QN_TYPE currqnumber = qnumbers[i];
    for (j = 0; j < 3; ++j, ++count)
    {
      indexes[count] = currqnumber % maxdims[count];
      currqnumber      = currqnumber / maxdims[count];
    }
    assert(currqnumber == 0);
  }
}

static int get_compressed_instructions(const int * const instructions, const int nr_instructions, 
    const int * const prev_hss, const int site, const int * const new_hss, int ** const compr_instr,
    int ** const compr_hss)
{
  int instr;
  int nrinstr = 0;
  *compr_instr = safe_malloc(nr_instructions * 3, int);
  *compr_hss = safe_malloc(nr_instructions * 3, int);

  /* Now loop over the different instructions and only execute the unique ones. */
  for (instr = 0; instr < nr_instructions; ++instr)
  {
    const int prevoperator = instructions[instr * 3 + 0];
    const int siteoperator = instructions[instr * 3 + 1];
    const int nextoperator = instructions[instr * 3 + 2];

    /* This instruction is the same as the previous one, thus already executed, just skip it. */
    /* If it is the first instruction, you have to execute it for sure */
    if (instr != 0 && prevoperator == instructions[(instr - 1) * 3 + 0] &&
        siteoperator == instructions[(instr - 1) * 3 + 1] && 
        new_hss[nextoperator] == new_hss[instructions[(instr - 1) * 3 + 2]])
      continue;

    /* save the new instructions */
    (*compr_instr)[3 * nrinstr + 0] = instructions[instr * 3 + 0];
    (*compr_instr)[3 * nrinstr + 1] = instructions[instr * 3 + 1];
    (*compr_instr)[3 * nrinstr + 2] = instructions[instr * 3 + 2];

    (*compr_hss)[3 * nrinstr + 0] = prev_hss[prevoperator];
    (*compr_hss)[3 * nrinstr + 1] = symsec_siteop(siteoperator, site);
    (*compr_hss)[3 * nrinstr + 2] = new_hss[nextoperator];
    ++nrinstr;
  }

  *compr_instr = realloc(*compr_instr, nrinstr * 3 * sizeof(int));
  *compr_hss   = realloc(*compr_hss, nrinstr * 3 * sizeof(int));
  return nrinstr;
}

static QN_TYPE change_new_to_old_ops(const QN_TYPE new, const int dimold, const int dimnew, int * 
    const oldtonew)
{
  /* [bra, ket, MPO] => [bra', ket', MPO] */
  assert(dimold >= dimnew);
  const int branew = new % dimnew;
  const QN_TYPE temp = new / dimnew;
  const int ketnew = temp % dimnew;
  const int MPO = temp / dimnew;
  int braold;
  int ketold;
  for (braold = branew; braold < dimold; ++braold)
    if (oldtonew[braold] == branew)
      break;
  assert(braold != dimold);
  for (ketold = ketnew; ketold < dimold; ++ketold)
    if (oldtonew[ketold] == ketnew)
      break;
  assert(ketold != dimold);

  return braold + ketold * dimold + MPO * dimold * dimold;
}

static QN_TYPE change_old_to_new_site(const QN_TYPE old, const int dimold, const int maxdims[],
    int * const oldtonew, const int is_left)
{
  if (!is_left)
  {
    assert(oldtonew[old % dimold] != -1);
    return oldtonew[old % dimold] + (old / dimold) * maxdims[0];
  }
  else
  {
    const QN_TYPE AB = maxdims[0] * maxdims[1];
    assert(oldtonew[old / AB] != -1);
    return old % AB + oldtonew[old / AB] * AB;
  }
}

static int * make_oldtonew(const struct symsecs * const internalss, const int bond)
{
  struct symsecs newss;
  int * result = safe_malloc(internalss->nrSecs, int);
  int cnt = 0;
  int i = 0;
  get_symsecs(&newss, bond);
  assert(internalss->nrSecs >= newss.nrSecs);

  for (i = 0; i < newss.nrSecs; ++i)
  {
    for (; cnt < internalss->nrSecs; ++cnt)
    {
      int j;
      for (j = 0; j < bookie.nrSyms; ++j)
        if (internalss->irreps[cnt][j] != newss.irreps[i][j])
          break;
      if (j == bookie.nrSyms)
      {
        result[cnt] = i;
        ++cnt;
        break;
      }
      else
        result[cnt] = -1;
    }
    assert(cnt < internalss->nrSecs || 
        (cnt == internalss->nrSecs && i == newss.nrSecs - 1));
  }
    assert(cnt < internalss->nrSecs || 
        (cnt == internalss->nrSecs && i == newss.nrSecs));
  for (; cnt < internalss->nrSecs; ++cnt)
    result[cnt] = -1;

  clean_symsecs(&newss, bond);
  return result;
}

#ifdef DEBUG
static int consistency_check_for_update_physical(const struct rOperators * const rops, 
    const struct siteTensor * const tens)
{
  /* I will check with this if the siteTensor is indeed a 1 site tensor and rops is a physical one. 
   * Also that the qnumberbonds and the indices of both tens and rops correspond. */
  int indicestens[3];
  int qnumbertens[3];
  int indicesrops[7];
  int qnumberrops[9];
  int i;
  if (tens->nrsites != 1)
  {
    fprintf(stderr, "%s@%s: nrsites of tens is not equal to 1 but %d.\n", __FILE__, __func__, 
        tens->nrsites);
    return 0;
  }
  if (rops->P_operator != 1)
  {
    fprintf(stderr, "%s@%s: rops is not a P_operator.\n", __FILE__, __func__);
    return 0;
  }

  siteTensor_give_indices(tens, indicestens);
  siteTensor_give_qnumberbonds(tens, qnumbertens);
  rOperators_give_indices(rops, indicesrops);
  rOperators_give_qnumberbonds(rops, qnumberrops);
  for (i = 0; i < 3; ++i)
  {
    if (!are_bra_and_ket_bonds(indicesrops[i], indicestens[i]) || 
        indicestens[i] != indicesrops[i + 3]) 
    {
      fprintf(stderr, "%s@%s: Something wrong with the indices array.\n", __FILE__, __func__);
      return 0;
    }
    if (!are_bra_and_ket_bonds(qnumberrops[i], qnumbertens[i]) || 
        qnumbertens[i] != qnumberrops[i + 3]) 
    {
      fprintf(stderr, "%s@%s: Something wrong with the qnumber array.\n", __FILE__, __func__);
      return 0;
    }
  }
  return 1;
}
#endif

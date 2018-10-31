#include <stdio.h>
#include <stdlib.h>

#include "instructions.h"
#include "instructions_qc.h"
#include "instructions_nn_hubbard.h"
#include "hamiltonian.h"
#include "sort.h"
#include "macros.h"
#include "network.h"
#include "debug.h"
#include "bookkeeper.h"

/* ========================================================================== */
/* ==================== DECLARATION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void sort_instructions(struct instructionset * const instructions);

static void sort_instructionsx(int ** instructions, double ** prefactors, 
                               const int nr_instructions, const int step);

static void print_DMRG_instructions(int * const instructions, 
                                    double * const prefactors, int * const hss, 
                                    const int nr_instructions, const int bond, 
                                    const int is_left);

static void print_T3NS_instructions(int * const instructions, 
                                    double * const prefactors, int * const hss, 
                                    const int nr_instructions, const int bond, 
                                    const int is_left);

static void print_merge_instructions(int * const instructions, 
                                     double * const prefactors, 
                                     const int nr_instructions, const int bond);

/* ========================================================================== */

void fetch_pUpdate(int ** const instructions, double ** const prefactors, 
                   int ** const hamsymsecs_of_new, int * const nr_instructions, 
                   const int bond, const int is_left)
{
        struct instructionset instr;
        switch(ham) {
        case QC :
                QC_fetch_pUpdate(&instr, bond, is_left);
                *instructions = instr.instr;
                *prefactors = instr.pref;
                *hamsymsecs_of_new = instr.hss_of_new;
                *nr_instructions = instr.nr_instr;
                break;
        case NN_HUBBARD :
                NN_H_fetch_pUpdate(instructions, prefactors, hamsymsecs_of_new, 
                                   nr_instructions, bond, is_left);
                break;
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n", 
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
        sort_instructionsx(instructions, prefactors, *nr_instructions, 3);
}

void fetch_bUpdate(struct instructionset * const instructions, const int bond, 
                   const int isleft)
{
        switch(ham) {
        case QC :
                QC_fetch_bUpdate(instructions, bond, isleft);
                break;
        case NN_HUBBARD :
                NN_H_fetch_bUpdate(instructions, bond, isleft);
                break;
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n", 
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
        sort_instructions(instructions);
}

void fetch_merge(int ** const instructions, int * const nr_instructions, 
                 double** const prefactors, const int bond)
{
        struct instructionset instr;
        switch(ham)
        {
        case QC :
                QC_fetch_merge(&instr, bond);
                *instructions = instr.instr;
                *prefactors = instr.pref;
                *nr_instructions = instr.nr_instr;
                break;
        case NN_HUBBARD :
                NN_H_fetch_merge(instructions, nr_instructions, 
                                 prefactors, bond);
                break;
        default:
                fprintf(stderr, "%s@%s: Unrecognized Hamiltonian.\n", 
                        __FILE__, __func__);
                exit(EXIT_FAILURE);
        }
}

void sortinstructions_toMPOcombos(int ** const instructions, 
                                  int ** const instrbegin, 
                                  double ** const prefactors, 
                                  const int nr_instructions, const int step, 
                                  int * const hss_of_Ops[step], 
                                  int ** const MPOinstr, int * const nrMPOinstr)
{
        int * temp = safe_malloc(nr_instructions, int); 
        int * newinstructions = safe_malloc(nr_instructions * step , int);
        double * newpref = safe_malloc(nr_instructions, double);
        int * idx;
        int i;
        const int hssdim = get_nr_hamsymsec();

        for (i = 0; i < nr_instructions; ++i)
        {
                int j;
                temp[i] = 0;
                for (j = step - 1; j >= 0; --j)
                        temp[i] = hss_of_Ops[j][(*instructions)[step * i + j]] + 
                                temp[i] * hssdim;
        }

        idx = quickSort(temp, nr_instructions);

        *instrbegin = safe_malloc(nr_instructions + 1, int);
        *MPOinstr   = safe_malloc(nr_instructions, int);
        *nrMPOinstr = 0;

        (*instrbegin)[(*nrMPOinstr)] = 0;
        (*MPOinstr)  [(*nrMPOinstr)] = temp[idx[0]];
        ++(*nrMPOinstr);
        newpref[0] = (*prefactors)[idx[0]];
        for (i = 0; i < step; ++i)
                newinstructions[0 * step + i] = 
                        (*instructions)[idx[0] * step + i];
        for (i = 1; i < nr_instructions; ++i) {
                int j;
                assert((*MPOinstr)[(*nrMPOinstr) - 1] <= temp[idx[i]]);

                newpref[i] = (*prefactors)[idx[i]];
                for (j = 0; j < step; ++j)
                        newinstructions[i * step + j] = 
                                (*instructions)[idx[i] * step + j];

                if ((*MPOinstr)[(*nrMPOinstr) - 1] != temp[idx[i]]) {
                        (*instrbegin)[(*nrMPOinstr)] = i;
                        (*MPOinstr)  [(*nrMPOinstr)] = temp[idx[i]];
                        ++(*nrMPOinstr);
                }
        }
        (*instrbegin)[(*nrMPOinstr)] = i;
        *instrbegin = realloc(*instrbegin, (*nrMPOinstr + 1) * sizeof(int));
        *MPOinstr   = realloc(*MPOinstr, *nrMPOinstr * sizeof(int));
        safe_free(*instructions);
        safe_free(*prefactors);
        safe_free(idx);
        safe_free(temp);
        *instructions = newinstructions;
        *prefactors = newpref;
}

void destroy_instructionset(struct instructionset * const instructions)
{
        safe_free(instructions->instr);
        safe_free(instructions->pref);
        safe_free(instructions->hss_of_new);
}

int get_next_unique_instr(int * const curr_instr, 
                          const struct instructionset * const instructions)
{
        /* instructions->instr is of the form:
         *   old1, old2, new
         * old1 and old2 should be different from prev instruction.
         * instructions->hss_of_new should also be different from prev
         * instruction. */
        if (*curr_instr == -1) {
                ++*curr_instr;
                return 1;
        } else {
                const int step = instructions->step;
                const int old1 = instructions->instr[step * *curr_instr];
                const int old2 = instructions->instr[step * *curr_instr + 1];
                const int old3 = instructions->instr[step * *curr_instr + 2];
                const int hss_old = instructions->hss_of_new[old3];

                for (++*curr_instr; *curr_instr < instructions->nr_instr; ++*curr_instr) {
                        const int new1 = instructions->instr[step * *curr_instr];
                        const int new2 = instructions->instr[step * *curr_instr + 1];
                        const int new3 = instructions->instr[step * *curr_instr + 2];
                        const int hss_new = instructions->hss_of_new[new3];
                        if (old1 != new1 || old2 != new2 || hss_old != hss_new)
                                return 1;
                }
                return 0;
        }
}

void print_instructions(int * const instructions, double * const prefactors, 
                        int * const hss, const int nr_instructions, 
                        const int bond, const int is_left, const char kind)
{
        switch(kind) {
        case 'd':
                print_DMRG_instructions(instructions, prefactors, hss, 
                                        nr_instructions, bond, is_left);
                break;
        case 'm':
                print_merge_instructions(instructions, prefactors, 
                                         nr_instructions, bond);
                break;
        case 't':
                print_T3NS_instructions(instructions, prefactors, hss, 
                                        nr_instructions, bond, is_left);
                break;
        default:
                fprintf(stderr, "%s@%s: Unknown option (%c)\n", __FILE__, __func__, kind);
        }
}

/* ========================================================================== */
/* ===================== DEFINITION STATIC FUNCTIONS ======================== */
/* ========================================================================== */

static void sort_instructions(struct instructionset * const instructions)
{
        const int step = instructions->step;
        const int nr_instr = instructions->nr_instr;
        int max[step];
        int *idx;
        int *instr_new  = safe_malloc(nr_instr * step, int);
        int *array      = safe_malloc(nr_instr, int);
        double *prefnew = instructions->pref == NULL ? 
                NULL : safe_malloc(nr_instr, double);
        int i, j;

        for (i = 0; i < step; ++i) {
                max[i] = -1;
                for (j = 0; j < nr_instr; ++j)
                        max[i] = max[i] < instructions->instr[j * step + i] + 1 ? 
                                instructions->instr[j * step + i] + 1 : max[i];
        }
        for (i = step - 2; i >= 0; --i) max[i] *= max[i + 1];
        for (i = 0; i < step - 1; ++i) max[i]  = max[i + 1];
        max[step - 1] = 1;

        for (i = 0; i < nr_instr; ++i) {
                array[i] = 0;
                for (j = 0; j < step; ++j)
                        array[i] += max[j] * instructions->instr[i * step + j];
        }

        idx = quickSort(array, nr_instr);
        for (i = 0; i < nr_instr; i++) {
                for (j = 0; j < step; ++j)
                        instr_new[i * step + j] = instructions->instr[idx[i] * step + j];

                if (prefnew != NULL) prefnew[i] = instructions->pref[idx[i]];
        }
        safe_free(instructions->instr);
        safe_free(instructions->pref);
        safe_free(array);
        safe_free(idx);
        instructions->instr      = instr_new;
        instructions->pref       = prefnew;
}

static void sort_instructionsx(int ** instructions, double ** prefactors, 
                               const int nr_instructions, const int step)
{
        int max[step];
        int *idx;
        int *instr_new  = safe_malloc(nr_instructions * step, int);
        int *array      = safe_malloc(nr_instructions, int);
        double *prefnew = prefactors == NULL ? NULL : safe_malloc(nr_instructions, double);
        int i, j;

        for (i = 0; i < step; ++i)
        {
                max[i] = -1;
                for (j = 0; j < nr_instructions; ++j)
                        max[i] = (max[i] < (*instructions)[j * step + i]+1) ? (*instructions)[j*step + i]+1 : max[i];
        }
        for (i = step - 2; i >= 0; --i) max[i] *= max[i + 1];
        for (i = 0; i < step - 1; ++i) max[i]  = max[i + 1];
        max[step - 1] = 1;

        for (i = 0; i < nr_instructions; ++i)
        {
                array[i] = 0;
                for (j = 0; j < step; ++j)
                        array[i] += max[j] * (*instructions)[i * step + j];
        }
        idx = quickSort(array, nr_instructions);
        for (i = 0; i < nr_instructions; i++)
        {
                for (j = 0; j < step; ++j)
                        instr_new[i * step + j] = (*instructions)[idx[i] * step + j];

                if (prefnew != NULL) prefnew[i] = (*prefactors)[idx[i]];
        }
        safe_free(*instructions);
        safe_free(*prefactors);
        safe_free(array);
        safe_free(idx);
        *instructions     = instr_new;
        *prefactors       = prefnew;
}

static void print_DMRG_instructions(int * const instructions, 
                                    double * const prefactors, 
                                    int * const hss, const int nr_instructions, 
                                    const int bond, const int is_left)
{
        const int site = netw.sitetoorb[netw.bonds[bond][is_left]];
        int bonds[3];
        int i;
        struct symsecs MPO;
        get_symsecs(&MPO, -1);

        get_bonds_of_site(netw.bonds[bond][is_left], bonds);
        assert(bond == bonds[2 * !is_left]);

        printf("================================================================================\n" 
               "Printing DMRG instructions for bond %d going %s.\n", bond, is_left ? "left" : "right");

        for (i = 0; i < nr_instructions; ++i)
        {
                char buffer[255];
                get_string_of_rops(buffer, instructions[i * 3 + 0], bond, is_left, 'e');
                printf("%14.8g * %-16s + ", prefactors[i], buffer);
                get_string_of_siteops(buffer, instructions[i * 3 + 1], site);
                printf("%-32s --> ", buffer);
                get_sectorstring(&MPO, hss[instructions[i * 3 + 2]], buffer);
                printf("(%s)", buffer);
                get_string_of_rops(buffer, instructions[i * 3 + 2], bonds[2 * is_left], is_left, 'c');
                printf("\t%s\n", buffer);
        }

        printf("================================================================================\n");
        clean_symsecs(&MPO, -1);
}

static void print_T3NS_instructions(int * const instructions, 
                                    double * const prefactors, int * const hss, 
                                    const int nr_instructions, const int bond, 
                                    const int is_left)
{
        int bonds[3];
        int i;
        int bond1, bond2, left1, left2;
        struct symsecs MPO;
        get_symsecs(&MPO, -1);

        assert(!is_psite(netw.bonds[bond][!is_left]));
        get_bonds_of_site(netw.bonds[bond][!is_left], bonds);
        assert((bonds[0] == bond && !is_left) || (bonds[1] == bond && !is_left) || 
               (bonds[2] == bond && is_left));

        if (bonds[0] == bond && !is_left) {
                bond1 = bonds[1];
                bond2 = bonds[2];
                left1 = 1;
                left2 = 0;
        } else if (bonds[1] == bond && !is_left) {
                bond1 = bonds[0];
                bond2 = bonds[2];
                left1 = 1;
                left2 = 0;
        } else {
                bond1 = bonds[0];
                bond2 = bonds[1];
                left1 = 1;
                left2 = 1;
        }

        printf("================================================================================\n" 
               "Printing T3NS update  instructions for bond %d going %s.\n", bond, 
               is_left ? "left" : "right");

        for (i = 0; i < nr_instructions; ++i) {
                char buffer[255];
                get_string_of_rops(buffer, instructions[i * 3 + 0], bond1, left1, 'e');
                printf("%14.8g * %-16s + ", prefactors[i], buffer);
                get_string_of_rops(buffer, instructions[i * 3 + 1], bond2, left2, 'e');
                printf("%-32s --> ", buffer);
                get_sectorstring(&MPO, hss[instructions[i * 3 + 2]], buffer);
                printf("(%s)", buffer);
                get_string_of_rops(buffer, instructions[i * 3 + 2], bond, is_left, 'c');
                printf("\t%s\n", buffer);
        }

        clean_symsecs(&MPO, -1);
}

static void print_merge_instructions(int * const instructions, 
                                     double * const prefactors, 
                                     const int nr_instructions, const int bond)
{
        const int isdmrg = is_dmrg_bond(bond);
        int i;
        const int step = 2 + !isdmrg;
        int bonds[step];
        int isleft[step];
        struct symsecs MPO;
        get_symsecs(&MPO, -1);

        if (isdmrg)
        {
                bonds[0] = bond;
                bonds[1] = bond;
                isleft[0] = 1;
                isleft[1] = 0;
        }
        else
        {
                int branching_site = netw.bonds[bond][is_psite(netw.bonds[bond][0])];
                assert(!is_psite(branching_site));
                get_bonds_of_site(branching_site, bonds);

                isleft[0] = 1;
                isleft[1] = 1;
                isleft[2] = 0;
        }

        printf("================================================================================\n" 
               "Printing merge instructions for bond %d.\n", bond);

        for (i = 0; i < nr_instructions; ++i)
        {
                char buffer[255];
                int j;
                printf("%14.8g * ", prefactors[i]);
                for (j = 0; j < step; ++ j)
                {
                        get_string_of_rops(buffer, instructions[i * step + j], bonds[j], isleft[j], 'e');
                        printf("%-16s%s", buffer, j == step - 1 ? "\n" : " + ");
                }
        }
        clean_symsecs(&MPO, -1);
}

#ifndef HAMILTONIAN_QC_H 
# define HAMILTONIAN_QC_H

#include "bookkeeper.h"

/**
 * \file hamiltonian_qc.h
 * \brief Implementation for the quantum chemistry hamiltonian with U1 SU2 and point-group 
 * symmetries.
 *
 * enum in hamiltonian wrapper is QC, QCSU2
 */

void QC_destroy_hamiltonian(void);

void QC_make_hamiltonian(char hamiltonianfile[], int su2);

void QC_get_physsymsecs(struct symsecs *res, int site);

void QC_get_hamiltoniansymsecs(struct symsecs * const res, const int bond);

int QC_get_nr_hamsymsec(void);

int QC_get_trivialhamsymsec(void);

int QC_hermitian_symsec(const int orig_symsec);

int QC_consistencynetworkinteraction(void);

double QC_el_siteop(const int siteop, const int braindex, const int ketindex);

double get_core(void);

void QC_tprods_ham(int * const nr_of_prods, int ** const possible_prods, 
                   const int resulting_symsec, const int site);

int QC_MPO_couples_to_singlet(const int n, const int MPO[n]);

void make_site_opType(int ** begin_opType, int **** tags_opType);

void QC_get_string_of_rops(char buffer[], const int ropsindex, const int bond, 
                           const int is_left, const char o);

void QC_get_string_of_siteops(char buffer[], const int siteindex, 
                              const int site);
#endif

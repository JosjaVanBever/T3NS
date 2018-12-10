#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <hdf5.h>
#include <unistd.h>
#include <omp.h>

#include "io_to_disk.h"
#include "sparseblocks.h"
#include "siteTensor.h"
#include "rOperators.h"
#include "bookkeeper.h"
#include "symsecs.h"
#include "network.h"
#include "macros.h"
#include "debug.h"
#include "hamiltonian.h"

static void write_symsec_to_disk(const hid_t id, const struct symsecs * const 
                                 ssec, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./symsec_%d", nmbr);
        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "nrSecs", &ssec->nrSecs, 1, THDF5_INT);
        write_attribute(group_id, "totaldims", &ssec->totaldims, 1, THDF5_INT);
        write_dataset(group_id, "./dims", ssec->dims, ssec->nrSecs, THDF5_INT);
        write_dataset(group_id, "./irreps", ssec->irreps, 
                      ssec->nrSecs * MAX_SYMMETRIES, THDF5_INT);
        write_dataset(group_id, "./fcidims", ssec->fcidims, ssec->nrSecs, 
                      THDF5_DOUBLE);
        H5Gclose(group_id);
}

static void read_symsec_from_disk(const hid_t id, struct symsecs * const ssec, 
                                  const int nmbr, int offset)
{
        char buffer[255];
        sprintf(buffer, "./symsec_%d", nmbr);
        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        read_attribute(group_id, "nrSecs", &ssec->nrSecs, THDF5_INT);
        read_attribute(group_id, "totaldims", &ssec->totaldims, THDF5_INT);

        ssec->dims = safe_malloc(ssec->nrSecs, int);
        read_dataset(group_id, "./dims", ssec->dims, THDF5_INT);

        int * dummyirreps = safe_malloc(ssec->nrSecs * offset, *dummyirreps);
        ssec->irreps = safe_malloc(ssec->nrSecs, *ssec->irreps);
        read_dataset(group_id, "./irreps", dummyirreps, THDF5_INT);

        for (int i = 0; i < ssec->nrSecs; ++i) {
                for (int j = 0; j < bookie.nrSyms; ++j) {
                        ssec->irreps[i][j] = dummyirreps[i * offset + j];
                }
        }
        safe_free(dummyirreps);
        
        ssec->fcidims = safe_malloc(ssec->nrSecs, double);
        read_dataset(group_id, "./fcidims", ssec->fcidims, THDF5_DOUBLE);
        H5Gclose(group_id);
}

static void write_bookkeeper_to_disk(const hid_t file_id)
{
        const hid_t group_id = H5Gcreate(file_id, "/bookkeeper", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        const int mxsyms = MAX_SYMMETRIES;
        write_attribute(group_id, "nrSyms", &bookie.nrSyms, 1, THDF5_INT);
        write_attribute(group_id, "Max_symmetries", &mxsyms, 1, THDF5_INT);
        write_attribute(group_id, "sgs", (int *) bookie.sgs,
                        bookie.nrSyms, THDF5_INT);
        write_attribute(group_id, "target_state", bookie.target_state, 
                        bookie.nrSyms, THDF5_INT);
        write_attribute(group_id, "nr_bonds", &bookie.nr_bonds, 1, THDF5_INT);

        for (int i = 0 ; i < bookie.nr_bonds; ++i) 
                write_symsec_to_disk(group_id, &bookie.list_of_symsecs[i], i);

        H5Gclose(group_id);
}

static void read_bookkeeper_from_disk(const hid_t file_id)
{
        const hid_t group_id = H5Gopen(file_id, "/bookkeeper", H5P_DEFAULT);

        read_attribute(group_id, "nrSyms", &bookie.nrSyms, THDF5_INT);
        if (bookie.nrSyms > MAX_SYMMETRIES) {
                fprintf(stderr, "Error: This wave function can not be read.\n"
                        "The currently compiled program can not run with the specified number of symmetries.\n"
                        "Specified: %d, maximal allowed: %d.\n"
                        "Recompile with a higher MAX_SYMMETRIES", bookie.nrSyms, MAX_SYMMETRIES);
                exit(EXIT_FAILURE);
        }
        int offset;
        read_attribute(group_id, "Max_symmetries", &offset, THDF5_INT);

        bookie.sgs = safe_malloc(bookie.nrSyms, enum symmetrygroup);
        read_attribute(group_id, "sgs", (int *) bookie.sgs, THDF5_INT);

        bookie.target_state = safe_malloc(bookie.nrSyms, int);
        read_attribute(group_id, "target_state", bookie.target_state, THDF5_INT);

        read_attribute(group_id, "nr_bonds", &bookie.nr_bonds, THDF5_INT);

        bookie.list_of_symsecs = safe_malloc(bookie.nr_bonds, struct symsecs);

        for (int i = 0 ; i < bookie.nr_bonds; ++i) 
                read_symsec_from_disk(group_id, &bookie.list_of_symsecs[i], i, 
                                      offset);

        H5Gclose(group_id);
}

static void write_sparseblocks_to_disk(const hid_t id, 
                                       const struct sparseblocks * block, 
                                       const int nrblocks, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./block_%d", nmbr);
        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "nrBlocks", &nrblocks, 1, THDF5_INT);
        write_dataset(group_id, "./beginblock", block->beginblock,
                      nrblocks + 1, THDF5_INT);
        write_dataset(group_id, "./tel", block->tel, 
                      block->beginblock[nrblocks], THDF5_EL_TYPE);

        H5Gclose(group_id);
}

static void read_sparseblocks_from_disk(const hid_t id, 
                                       struct sparseblocks * block, 
                                       const int nrblocks, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./block_%d", nmbr);
        int bloccount;

        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        read_attribute(group_id, "nrBlocks", &bloccount, THDF5_INT);
        assert(bloccount == nrblocks);

        block->beginblock = safe_malloc(nrblocks + 1, int);
        read_dataset(group_id, "./beginblock", block->beginblock, THDF5_INT);

        block->tel = safe_malloc(block->beginblock[nrblocks], EL_TYPE);
        read_dataset(group_id, "./tel", block->tel, THDF5_EL_TYPE);

        H5Gclose(group_id);
}

static void write_siteTensor_to_disk(const hid_t id, const struct siteTensor * 
                                     const tens, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./tensor_%d", nmbr);
        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "nrsites", &tens->nrsites, 1, THDF5_INT);
        write_attribute(group_id, "sites", tens->sites, tens->nrsites, THDF5_INT);
        write_attribute(group_id, "nrblocks", &tens->nrblocks, 1, THDF5_INT);
        write_dataset(group_id, "./qnumbers", tens->qnumbers, 
                      tens->nrblocks * tens->nrsites, THDF5_QN_TYPE);
        write_sparseblocks_to_disk(group_id, &tens->blocks, tens->nrblocks, 0);
        H5Gclose(group_id);
}

static void read_siteTensor_from_disk(const hid_t id, struct siteTensor * 
                                     const tens, const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./tensor_%d", nmbr);
        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        read_attribute(group_id, "nrsites", &tens->nrsites, THDF5_INT);

        tens->sites = safe_malloc(tens->nrsites, int);
        read_attribute(group_id, "sites", tens->sites, THDF5_INT);
        read_attribute(group_id, "nrblocks", &tens->nrblocks, THDF5_INT);

        tens->qnumbers = safe_malloc(tens->nrblocks * tens->nrsites, QN_TYPE);
        read_dataset(group_id, "./qnumbers", tens->qnumbers, THDF5_QN_TYPE);
        read_sparseblocks_from_disk(group_id, &tens->blocks, tens->nrblocks, 0);
        H5Gclose(group_id);
}

static void write_T3NS_to_disk(const hid_t file_id, 
                               const struct siteTensor * const T3NS)
{
        const hid_t group_id = H5Gcreate(file_id, "/T3NS", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "nrSites", &netw.sites, 1, THDF5_INT);

        for(int i = 0 ; i < netw.sites; ++i) 
                write_siteTensor_to_disk(group_id, &T3NS[i], i);

        H5Gclose(group_id);
}

static void read_T3NS_from_disk(const hid_t file_id, 
                                struct siteTensor ** const T3NS)
{
        const hid_t group_id = H5Gopen(file_id, "/T3NS", H5P_DEFAULT);
        int nrsit;

        read_attribute(group_id, "nrSites", &nrsit, THDF5_INT);
        assert(nrsit == netw.sites);

        *T3NS = safe_malloc(nrsit, struct siteTensor);

        for(int i = 0 ; i < netw.sites; ++i) 
                read_siteTensor_from_disk(group_id, &(*T3NS)[i], i);

        H5Gclose(group_id);
}

static void write_rOperator_to_disk(const hid_t id,
                                    const struct rOperators * const rOp,
                                    const int nmbr)
{
        char buffer[255];
        if (rOp->is_left == -1 || rOp->bond_of_operator == -1)
                return;

        sprintf(buffer, "./rOperator_%d", nmbr);

        hid_t group_id = H5Gcreate(id, buffer, H5P_DEFAULT, 
                                   H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "bond_of_operator", &rOp->bond_of_operator, 1, THDF5_INT);
        write_attribute(group_id, "is_left", &rOp->is_left, 1, THDF5_INT);
        write_attribute(group_id, "P_operator", &rOp->P_operator, 1, THDF5_INT);
        write_attribute(group_id, "nrhss", &rOp->nrhss, 1, THDF5_INT);

        write_dataset(group_id, "./begin_blocks_of_hss",
                      rOp->begin_blocks_of_hss, rOp->nrhss + 1, THDF5_INT);
        write_dataset(group_id, "./qnumbers", rOp->qnumbers, 
                      rOp->begin_blocks_of_hss[rOp->nrhss] *
                      rOperators_give_nr_of_couplings(rOp), THDF5_QN_TYPE);

        write_attribute(group_id, "nrops", &rOp->nrops, 1, THDF5_INT);
        write_dataset(group_id, "./hss_of_ops", rOp->hss_of_ops, 
                      rOp->nrops, THDF5_INT);

        for (int i = 0; i < rOp->nrops; ++i) {
                const int nr_blocks = rOperators_give_nr_blocks_for_operator(rOp, i);
                write_sparseblocks_to_disk(group_id, &rOp->operators[i], nr_blocks, i);
        }
        H5Gclose(group_id);
}

static void read_rOperator_from_disk(const hid_t id,
                                    struct rOperators * const rOp, 
                                    const int nmbr)
{
        char buffer[255];
        sprintf(buffer, "./rOperator_%d", nmbr);

        const hid_t group_id = H5Gopen(id, buffer, H5P_DEFAULT);

        read_attribute(group_id, "bond_of_operator", &rOp->bond_of_operator, THDF5_INT);
        read_attribute(group_id, "is_left", &rOp->is_left, THDF5_INT);
        read_attribute(group_id, "P_operator", &rOp->P_operator, THDF5_INT);
        read_attribute(group_id, "nrhss", &rOp->nrhss, THDF5_INT);

        rOp->begin_blocks_of_hss = safe_malloc(rOp->nrhss + 1, int);
        read_dataset(group_id, "./begin_blocks_of_hss",
                     rOp->begin_blocks_of_hss, THDF5_INT);

        rOp->qnumbers = safe_malloc(rOp->begin_blocks_of_hss[rOp->nrhss] *
                              rOperators_give_nr_of_couplings(rOp), QN_TYPE);
        read_dataset(group_id, "./qnumbers", rOp->qnumbers, THDF5_QN_TYPE);

        read_attribute(group_id, "nrops", &rOp->nrops, THDF5_INT);

        rOp->hss_of_ops = safe_malloc(rOp->nrops, int);
        read_dataset(group_id, "./hss_of_ops", rOp->hss_of_ops, THDF5_INT);

        rOp->operators = safe_malloc(rOp->nrops, struct sparseblocks);
        for (int i = 0; i < rOp->nrops; ++i) {
                const int nr_blocks = rOperators_give_nr_blocks_for_operator(rOp, i);
                if (nr_blocks == 0) {
                        rOp->operators[i].beginblock = NULL;
                        rOp->operators[i].tel = NULL;
                }
                else
                        read_sparseblocks_from_disk(group_id, &rOp->operators[i], 
                                                    nr_blocks, i);
        }
        H5Gclose(group_id);
}

static void write_rOps_to_disk(const hid_t id,
                               const struct rOperators * const rOps)
{
        const hid_t group_id = H5Gcreate(id, "/rOps", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "nrOps", &netw.nr_bonds, 1, THDF5_INT);

        for (int i = 0 ; i < netw.nr_bonds; ++i) {
                write_rOperator_to_disk(group_id, &rOps[i], i);
        }

        H5Gclose(group_id);
}

static void read_rOps_from_disk(const hid_t id,
                               struct rOperators ** const rOps)
{
        const hid_t group_id = H5Gopen(id, "/rOps", H5P_DEFAULT);
        int nrbonds;

        read_attribute(group_id, "nrOps", &nrbonds, THDF5_INT);
        assert(nrbonds == netw.nr_bonds);

        *rOps = safe_malloc(nrbonds, struct rOperators);
        for (int i = 0 ; i < netw.nr_bonds; ++i) {
                read_rOperator_from_disk(group_id, &(*rOps)[i], i);
        }

        H5Gclose(group_id);
}

static void write_network_to_disk(const hid_t id)
{
        const hid_t group_id = H5Gcreate(id, "/network", H5P_DEFAULT, 
                                         H5P_DEFAULT, H5P_DEFAULT);

        write_attribute(group_id, "nr_bonds", &netw.nr_bonds, 1, THDF5_INT);
        write_dataset(group_id, "./bonds", netw.bonds, netw.nr_bonds * 2, THDF5_INT);

        write_attribute(group_id, "psites", &netw.psites, 1, THDF5_INT);
        write_attribute(group_id, "sites", &netw.sites, 1, THDF5_INT);
        write_dataset(group_id, "./sitetoorb", netw.sitetoorb, 
                      netw.sites, THDF5_INT);

        write_attribute(group_id, "sweeplength", &netw.sweeplength, 1, THDF5_INT);
        write_dataset(group_id, "./sweep", netw.sweep, 
                      netw.sweeplength, THDF5_INT);

        H5Gclose(group_id);
}

static void read_network_from_disk(const hid_t id)
{
        const hid_t group_id = H5Gopen(id, "/network", H5P_DEFAULT);

        read_attribute(group_id, "nr_bonds", &netw.nr_bonds, THDF5_INT);

        netw.bonds = malloc(netw.nr_bonds * sizeof netw.bonds[0]);
        read_dataset(group_id, "./bonds", (int*) netw.bonds, THDF5_INT);

        read_attribute(group_id, "psites", &netw.psites, THDF5_INT);
        read_attribute(group_id, "sites", &netw.sites, THDF5_INT);

        netw.sitetoorb = safe_malloc(netw.sites, int);
        read_dataset(group_id, "./sitetoorb", netw.sitetoorb, THDF5_INT);

        read_attribute(group_id, "sweeplength", &netw.sweeplength, THDF5_INT);

        netw.sweep = safe_malloc(netw.sweeplength, int);
        read_dataset(group_id, "./sweep", netw.sweep, THDF5_INT);

        H5Gclose(group_id);
        print_network();

        create_nr_left_psites();
        create_order_psites();
}

/* ========================================================================== */

static void make_h5f_name(const char * hdf5_loc, const char hdf5file[], 
                          const int size, char hdf5_resulting[size])
{
        strncpy(hdf5_resulting, hdf5_loc, size - 1);
        hdf5_resulting[size - 1] = '\0';

        int len = strlen(hdf5_resulting);
        if (hdf5_resulting[len - 1] != '/' && hdf5file[0] != '/') {
                hdf5_resulting[len] = '/';
                hdf5_resulting[len + 1] = '\0';
                ++len;
        }
        if (hdf5_resulting[len - 1] == '/' && hdf5file[0] == '/') {
                hdf5_resulting[len - 1] = '\0';
                --len;
        }

        strncat(hdf5_resulting, hdf5file, size - len);
        hdf5_resulting[size - 1] = '\0';
}

void write_to_disk(const char * hdf5_loc, const struct siteTensor * const T3NS, 
                   const struct rOperators * const ops)
{
        if (hdf5_loc == NULL)
                return;

        hid_t file_id;
        const char hdf5nam[] = "T3NScalc.h5";
        char hdf5file[MY_STRING_LEN];
        make_h5f_name(hdf5_loc, hdf5nam, MY_STRING_LEN, hdf5file);

        file_id = H5Fcreate(hdf5file, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        write_network_to_disk(file_id);
        write_bookkeeper_to_disk(file_id);
        write_hamiltonian_to_disk(file_id);
        write_T3NS_to_disk(file_id, T3NS);
        write_rOps_to_disk(file_id, ops);

        H5Fclose(file_id);
}

void read_from_disk(const char filename[], struct siteTensor ** const T3NS, 
                    struct rOperators ** const ops)
{
        if (access(filename, F_OK) != 0) {
                fprintf(stderr, "Error in %s: Can not read from disk.\n"
                        "%s was not found.\n", __func__, filename);
                exit(EXIT_FAILURE);
        }

        hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

        read_network_from_disk(file_id);
        read_bookkeeper_from_disk(file_id);
        read_hamiltonian_from_disk(file_id);
        read_T3NS_from_disk(file_id, T3NS);
        read_rOps_from_disk(file_id, ops);

        H5Fclose(file_id);
}

void write_attribute(hid_t group_id, const char atrname[], const void * atr, 
                     hsize_t size, enum hdf5type kind)
{
        hid_t datatype_arr[] = {
                H5T_STD_I32LE, H5T_IEEE_F64LE, EL_TYPE_H5, QN_TYPE_H5
        };

        if (atr == NULL)
                return;
        hid_t datatype = datatype_arr[kind];

        const hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        const hid_t attribute_id = H5Acreate(group_id, atrname, datatype, 
                                       dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

        H5Awrite(attribute_id, datatype, atr);

        H5Aclose(attribute_id);
        H5Sclose(dataspace_id);
}

void read_attribute(hid_t id, const char atrname[], void * atr, 
                    enum hdf5type kind)
{
        hid_t datatype_arr[] = {
                H5T_STD_I32LE, H5T_IEEE_F64LE, EL_TYPE_H5, QN_TYPE_H5
        };
        hid_t datatype = datatype_arr[kind];

        hid_t attribute_id = H5Aopen(id, atrname, H5P_DEFAULT);
        H5Aread(attribute_id, datatype, atr);
        H5Aclose(attribute_id);
}

void write_dataset(hid_t id, const char datname[], const void * dat, hsize_t size,
                   enum hdf5type kind)
{
        hid_t datatype_arr[] = {
                H5T_STD_I32LE, H5T_IEEE_F64LE, EL_TYPE_H5, QN_TYPE_H5
        };

        if (dat == NULL) { return; }
        hid_t datatype = datatype_arr[kind];

        hid_t dataspace_id = H5Screate_simple(1, &size, NULL);
        hid_t dataset_id = H5Dcreate(id, datname, datatype, dataspace_id, 
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        H5Dwrite (dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
}

void read_dataset(hid_t id, const char datname[], void * dat,
                  enum hdf5type kind)
{
        hid_t datatype_arr[] = {
                H5T_STD_I32LE, H5T_IEEE_F64LE, EL_TYPE_H5, QN_TYPE_H5
        };

        hid_t datatype = datatype_arr[kind];
        hid_t dataset_id = H5Dopen(id, datname, H5P_DEFAULT);
        H5Dread(dataset_id, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, dat);
        H5Dclose(dataset_id);
}

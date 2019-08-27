#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "mpi.h"

#define TOTAL_RUNS 100

#define MAX_RANDOM_NUM (1<<20)
#define MASTER 0
#define EPSILON 1e-9

#include "mmio-wrapper.h"
#include "util.h"
#include "partition.h"

MPI_Datatype procs_info_type;
enum tag {
    REQUEST_TAG, RECEIVE_TAG
};


/**
 * Creates two MPI derived datatypes for the respective structs
 */
void create_mpi_datatypes(MPI_Datatype *procs_info_type) {
    MPI_Datatype oldtypes[2];
    MPI_Aint offsets[2];
    int blockcounts[2];

    /* create `proc_info_t` datatype */
    offsets[0] = 0;
    oldtypes[0] = MPI_INT;
    blockcounts[0] = 9;

    MPI_Type_create_struct(1, blockcounts, offsets, oldtypes, procs_info_type);
    MPI_Type_commit(procs_info_type);
}


/**
 *
 * @param argc
 * @param argv
 * @return
 * Main method of the SpMV model
 * This model will able to predict the run time of the sparse matrix
 */
int main(int argc, char *argv[]) {

    char *in_file;
    int nRanks, *row_ptr, *col_ptr;
    double *val_ptr;


    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        in_file = argv[1];
        nRanks = argv[2];
    }
    int sqrRank = sqrt(nRanks);
    for (int rank = 0; rank < nRanks; ++rank) {
        printf("[%d] Start reading\n",rank);
        if (test_csr_read_2D_partitioned_mat(in_file, &row_ptr, &col_ptr, &val_ptr, sqrRank, rank) != 0) {
            fprintf(stderr, "read_matrix: failed\n");
            exit(EXIT_FAILURE);
        }
        printf("[%d] Done reading\n",rank);
    }
    /* MPI: end */

    return 0;
}

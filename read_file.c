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

double *matMull(int *row_ptr, int *col_ptr, double *val_ptr, double *x, int nRow, int firstRow, double *y) {

    /// multiplication
    for (int i = 0; i <nRow; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * x[col_ptr[k] - firstRow];
    }
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
        nRanks = atoi(argv[2]);
    }
    int sqrRank = sqrt(nRanks);
    printf("[%d] allocate space for vector\n", nRanks);
    double *y = (double *) calloc_or_exit(123761, sizeof(double));
    double *x = (double *) malloc_or_exit(123761 * sizeof(double));
    printf("[%d] Done allocating\n", nRanks);
    for (int rank = 0; rank < nRanks; ++rank) {
        char *file;
        stpcpy(in_file, file);
        printf("[%d] Start reading\n",rank);
        if (test_csr_read_2D_partitioned_mat(file, &row_ptr, &col_ptr, &val_ptr, sqrRank, rank) != 0) {
            fprintf(stderr, "read_matrix: failed\n");
            exit(EXIT_FAILURE);
        }
        printf("[%d] Done reading\n",rank);
        matMull(row_ptr, col_ptr, val_ptr, x, 123761, (rank/sqrRank)*123761, y);
        printf("[%d] Done Multiplication\n",rank);
    }
    /* MPI: end */

    return 0;
}


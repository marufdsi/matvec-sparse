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

double *matMull(int rank, int *row_ptr, int *col_ptr, double *val_ptr, double *x, int nRow, int firstRow, double *y) {

    /// multiplication
    for (int i = 0; i <nRow; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * x[col_ptr[k] - firstRow];
    }
}

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
    double comp_time = 0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    int total_run = 100, nRanks, rank;
    int *row_ptr, *col_ptr;
    double *val_ptr, *x, *y;
    proc_info_t *ranks_info;
    proc_info_t *procs_info;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    /// Create custom mpi data type
    create_mpi_datatypes(&procs_info_type);
    /// Initialize process info to share among the ranks
    procs_info = (proc_info_t *) malloc_or_exit(nRanks * sizeof(proc_info_t));
    ranks_info = (proc_info_t *) malloc_or_exit(nRanks * sizeof(proc_info_t));

    int sqrRank = sqrtl(nRanks);
    int col_rank = rank/sqrRank; //which col of proc am I
    int row_rank = rank%sqrRank; //which row of proc am I

    //initialize communicators
    MPI_Comm commrow;
    MPI_Comm_split (MPI_COMM_WORLD, row_rank, rank, &commrow);

    MPI_Comm commcol;
    MPI_Comm_split (MPI_COMM_WORLD, col_rank, rank, &commcol);

    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        in_file = argv[1];
    }
    if (csr_read_2D_partitioned_mat(in_file, &row_ptr, &col_ptr, &val_ptr, &ranks_info, sqrRank, rank) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }

    printf("[%d] Reading Matrix done!\n", rank);
    y = (double *) calloc_or_exit(ranks_info[rank].M, sizeof(double));
    x = (double *) malloc_or_exit(ranks_info[rank].M * sizeof(double));
    for (int i = 0; i < ranks_info[rank].M; ++i) {
        x[i] = 1.0;
    }

    printf("[%d] Vector creation done!\n", rank);

    /// Share process info among all the processes
    MPI_Allgather(&ranks_info[rank], 1, procs_info_type, procs_info, 1, procs_info_type, MPI_COMM_WORLD);

    MPI_Bcast (x, ranks_info[rank].M, MPI_FLOAT, col_rank, commcol);

    printf("[%d] Reading Matrix and broad casting done!\n", rank);
    /// Start sparse matrix vector multiplication for each rank
    double start_time = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    for (int r = 0; r < total_run; ++r) {
        //broadcast X along column communicator
        MPI_Bcast (x, ranks_info[rank].M, MPI_FLOAT, col_rank, commcol); //col_rank is the one with the correct information

        // Multiplication
        matMull(rank, row_ptr, col_ptr, val_ptr, x, ranks_info[rank].M, ranks_info[rank].first_row, y);

        //reduce Y along row communicator
        MPI_Reduce(y, x, ranks_info[rank].M, MPI_FLOAT, MPI_SUM, row_rank, commrow);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    comp_time = (MPI_Wtime() - start_time) * 1000.00;
    avg_time = comp_time / total_run;

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nRanks;

    int max_nnz = 0, max_row = 0;
    double max_nnz_per_row = 0.0, nnz_per_row = procs_info[rank].NZ/procs_info[rank].M;
    MPI_Reduce(&procs_info[rank].NZ, &max_nnz, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&procs_info[rank].M, &max_row, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&nnz_per_row, &max_nnz_per_row, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);

    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n",
               rank, min_time, max_time, mean, procs_info[rank].NZ);
        char *_ptr = strtok(in_file, "_");
        FILE *resultCSV;
        FILE *checkFile;
        if ((checkFile = fopen("CSR_Random_BrCast_Reduce_SpMV.csv", "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen("CSR_Random_BrCast_Reduce_SpMV.csv", "a"))) {
                fprintf(stderr, "fopen: failed to open file CSR_Random_BrCast_Reduce_SpMV.csv");
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen("CSR_Random_BrCast_Reduce_SpMV.csv", "w"))) {
                fprintf(stderr, "fopen: failed to open file CSR_Random_BrCast_Reduce_SpMV.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "Name,MatrixSize,MaxRow,MinTime,MaxTime,AvgTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock\n");
        }

        fprintf(resultCSV, "%s,%d,%d,%10.3lf,%10.3lf,%10.3lf,%d,%d,%10.3lf,%d\n", _ptr,procs_info[rank].N, max_row, min_time, max_time,
                mean, total_run, nRanks, max_nnz_per_row, max_nnz);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file CSR_Random_BrCast_Reduce_SpMV.csv");
            exit(EXIT_FAILURE);
        }
    }
    free(x);
    free(y);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


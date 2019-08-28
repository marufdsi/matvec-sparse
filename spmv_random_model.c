#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

#include "mpi.h"

#define TOTAL_RUNS 100

#define MAX_RANDOM_NUM (1<<20)
#define MASTER 0
#define EPSILON 1e-9

#include "mmio-wrapper.h"
#include "util.h"
#include "partition.h"


double *matMull(int rank, int *row_ptr, int *col_ptr, double *val_ptr, double *x, int nRow, int startCol, double *y) {

    /// multiplication
    for (int i = 0; i < nRow; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * x[col_ptr[k] - startCol];
    }
}

int createCSRMat(int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int nnz_per_block, int startCol, int rank){
    row_ptr = (int *) calloc_or_exit(mat_row+1, sizeof(int));
    col_ptr = (int *) malloc_or_exit(nnz_per_block * sizeof(int));
    val_ptr = (double *) malloc_or_exit(nnz_per_block * sizeof(double));
    srand(time(NULL) * (rank + 1));
    for (int k = 0; k < nnz_per_block; ++k) {
        int randColIdx = rand() % mat_row;
        int rowIdx = rank%2 == 0 ? (mat_row-1)-(k%mat_row) : (k%mat_row);
        row_ptr[rowIdx]++;
        col_ptr[k] = startCol + randColIdx;
        val_ptr[k] = 1.0;
    }
    for (int i = 0, cumsum = 0; i < mat_row; i++) {
        int temp = row_ptr[i];
        row_ptr[i] = cumsum;
        cumsum += temp;
    }
    row_ptr[mat_row] = nnz_per_block;
    return 0;
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
    double comp_time = 0, bcast_time = 0.0, matmul_time = 0.0, reduce_time = 0.0, min_time = 0.0, max_time = 0.0,
            avg_time = 0.0, mean = 0.0, avg_bcast_time = 0.0, avg_matmul_time = 0.0, avg_reduce_time = 0.0,
            *val_ptr, *x, *y;
    int total_run = 100, nRanks, rank, *row_ptr, *col_ptr, _size, mat_row, nnz_per_block;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int sqrRank = sqrt(nRanks);
    int row_rank = rank / sqrRank; //which col of proc am I
    int col_rank = rank % sqrRank; //which row of proc am I

    //initialize communicators
    MPI_Comm commrow;
    MPI_Comm_split(MPI_COMM_WORLD, row_rank, rank, &commrow);

    MPI_Comm commcol;
    MPI_Comm_split(MPI_COMM_WORLD, col_rank, rank, &commcol);

    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        _size = atoi(argv[1]);
        mat_row = atoi(argv[2]);
        nnz_per_block = atoi(argv[3]);
        if (argc > 4)
            total_run = atoi(argv[4]);
    }
    if (createCSRMat(row_ptr, col_ptr, val_ptr, mat_row, nnz_per_block, col_rank * mat_row, rank) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }

    y = (double *) calloc_or_exit(mat_row, sizeof(double));
    x = (double *) malloc_or_exit(mat_row * sizeof(double));
    for (int i = 0; i < mat_row; ++i) {
        x[i] = 1.0;
    }
    /// Start sparse matrix vector multiplication for each rank
    double start_bcast_time = 0.0, start_matmul_time = 0.0, start_reduce_time = 0.0;
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    for (int r = 0; r < total_run; ++r) {
        start_bcast_time = MPI_Wtime();
        //broadcast X along column communicator
        MPI_Bcast(x, mat_row, MPI_FLOAT, col_rank,
                  commcol); //col_rank is the one with the correct information
        bcast_time += (MPI_Wtime() - start_bcast_time) * 1000.00;

        start_matmul_time = MPI_Wtime();
        // Multiplication
        matMull(rank, row_ptr, col_ptr, val_ptr, x, mat_row, col_rank * mat_row, y);
        matmul_time += (MPI_Wtime() - start_matmul_time) * 1000.00;

        start_reduce_time = MPI_Wtime();
        //reduce Y along row communicator
        MPI_Reduce(y, x, mat_row, MPI_FLOAT, MPI_SUM, row_rank, commrow);
        reduce_time += (MPI_Wtime() - start_reduce_time) * 1000.00;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    comp_time = (MPI_Wtime() - start_time) * 1000.00;
    avg_time = comp_time / total_run;
    avg_bcast_time = bcast_time / total_run;
    avg_matmul_time = matmul_time / total_run;
    avg_reduce_time = reduce_time / total_run;

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nRanks;

    double nnz_per_row = (double)nnz_per_block/mat_row;
    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n",
               rank, min_time, max_time, mean, nnz_per_block);
        char *_ptr = strtok(in_file, "_");
        FILE *resultCSV;
        FILE *checkFile;
        if ((checkFile = fopen("CSR_SpMV_Model_of_Random_BrCast_Reduce.csv", "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen("CSR_SpMV_Model_of_Random_BrCast_Reduce.csv", "a"))) {
                fprintf(stderr, "fopen: failed to open file CSR_SpMV_Model_of_Random_BrCast_Reduce.csv");
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen("CSR_SpMV_Model_of_Random_BrCast_Reduce.csv", "w"))) {
                fprintf(stderr, "fopen: failed to open file CSR_SpMV_Model_of_Random_BrCast_Reduce.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "MatrixSize,PartitionRow,MinTime,MaxTime,AvgTime,AvgBcastTime,AvgMatmulTime,AvgReduceTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock\n");
        }

        fprintf(resultCSV, "%d,%d,%10.3lf,%10.3lf,%10.3lf,%10.3lf,%10.3lf,%10.3lf,%d,%d,%10.3lf,%d\n", _size, mat_row,
                min_time, max_time, mean, avg_bcast_time, avg_matmul_time, avg_reduce_time, total_run, nRanks, nnz_per_row, nnz_per_block);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file CSR_SpMV_Model_of_Random_BrCast_Reduce.csv");
            exit(EXIT_FAILURE);
        }
    }
    free(x);
    free(y);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


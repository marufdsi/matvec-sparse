#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include "mpi.h"

#undef DEBUG

#define TOTAL_RUNS 100

#define MAX_RANDOM_NUM (1<<20)
#define MASTER 0
#define EPSILON 1e-9

#include "mmio-wrapper.h"
#include "util.h"
#include "partition.h"


double *matMullComputationOnly(int rank, int *row_ptr, int *col_ptr, double *val_ptr, double *buf_x, int mat_row) {
    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(mat_row, sizeof(double));
    int first_row = rank*mat_row;
    /// Sparse Matrix Vector Multiplication without Communication
    for (int i = 0; i < mat_row; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            y[i] += val_ptr[k] * buf_x[col_ptr[k] - first_row];
        }
    }
    return y;
}

int main(int argc, char *argv[]) {

    double t, comp_time;
    int nprocs,     /* number of tasks/processes */
            rank;       /* id of task/process */

    int *row_ptr,     /* row index for all matrix elements */
            *col_ptr;     /* column index for all matrix elements */
    double *val_ptr, /* value for all matrix elements */
            *buf_x;      /* value for all x vector elements */

    /*******************************************/

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int mat_size = 0, nonZero = 0, nonZeroPerRow = 0, total_run = 100, mat_row = 0, bandwidth=0, gap = 1;

    if (argc < 4) {
        printf("Usage: %s matrix_row nonZeroPerRow [Total_Runs]\n", argv[0]);
        return 0;
    } else {
        mat_row = atoi(argv[1]);
        mat_size = atoi(argv[2]);
        nonZeroPerRow = atoi(argv[3]);
        if (argc > 4)
            gap = atoi(argv[4]);
        if (argc > 5)
            total_run = atoi(argv[5]);
    }

//    mat_size = mat_row*nprocs;
    nonZero = nonZeroPerRow*mat_row;
    if (nonZeroPerRow<=0){
        printf("[%d], There will must one non zero column in the matrix in every row\n", rank);
        return 0;
    }
    if(nonZeroPerRow > mat_row) {
        if(rank == MASTER) {
            printf("[%d] nonzero=%d, max nonzero=%d, number process=%d\n", rank, nonZeroPerRow, mat_row, nprocs);
        }
        nonZeroPerRow = mat_row;
        nonZero = nonZeroPerRow*mat_row;
    }
    if(nonZero <= 0){
        printf("[%d] Matrix can not be sized zero=%d\n", rank, nonZero);
    }

    row_ptr = (int *)malloc( (mat_row+1) * sizeof(int) );
    col_ptr = (int *)malloc( nonZero * sizeof(int) );
    val_ptr = (double *)malloc( nonZero * sizeof(double));
    /*if(csr_diagonal_mat(rank, row_ptr, col_ptr, val_ptr, mat_row, nonZeroPerRow) != 1){
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nprocs, mat_size, nonZeroPerRow);
    }*/

    printf("[%d] Gap=%d\n", rank, gap);
    if(csr_diagonal_mat_with_bandwidth(rank, row_ptr, col_ptr, val_ptr, mat_row, nonZeroPerRow, gap, &bandwidth) != 1){
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nprocs, mat_size, nonZeroPerRow);
    }

    printf("[%d] Bandwidth=%d\n", rank, bandwidth);
    buf_x = (double *) malloc_or_exit(mat_row * sizeof(double));
    for (int i = 0; i < mat_row; i++) {
        buf_x[i] = 1.00;
    }
    /* Matrix-vector multiplication for each processes */
    double totalTime = 0.0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    double *res;
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    for (int r = 0; r < total_run; ++r) {
        res = matMullComputationOnly(rank, row_ptr, col_ptr, val_ptr, buf_x, mat_row);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    totalTime = (MPI_Wtime() - t) * 1000.00;
    avg_time = totalTime / total_run;
//    printf("[%d] Matrix Multiplication done\n", rank);

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nprocs;

    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n",
               rank, min_time, max_time, mean, nonZero);
        FILE *resultCSV;
        FILE *checkFile;
        /// File Name CSR_MPI_SpMV.csv
        if ((checkFile = fopen("MPI_SpMV_Model.csv", "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen("MPI_SpMV_Model.csv", "a"))) {
                fprintf(stderr, "fopen: failed to open file MPI_SpMV_Model.csv");
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen("MPI_SpMV_Model.csv", "w"))) {
                fprintf(stderr, "fopen: failed to open file MPI_SpMV_Model.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "MatrixSize,MinTime,MaxTime,AvgTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock,Bandwidth,Gap\n");
        }

        fprintf(resultCSV, "%d,%10.3lf,%10.3lf,%10.3lf,%d,%d,%d,%d,%d,%d\n", mat_size, min_time, max_time, mean,
                total_run, nprocs, nonZeroPerRow, nonZero, bandwidth, gap);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file MPI_SpMV_Model");
            exit(EXIT_FAILURE);
        }
    }

    free(row_ptr);
    free(col_ptr);
    free(val_ptr);
    free(buf_x);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


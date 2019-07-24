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


double *matMullComputationOnly(int rank, int *buf_i_idx, int *buf_j_idx, double *buf_values, double *buf_x, int first_row, int mat_row, int nonZero) {
    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(mat_row, sizeof(double));
    /// Sparse Matrix Vector Multiplication without Communication
    for (int k = 0; k < nonZero; k++) {
//        if((buf_i_idx[k] - first_row) >= mat_row){
//            printf("[%d] row Outof index for i=%d, j=%d, val=%lf, first row=%d\n", rank, buf_i_idx[k], buf_j_idx[k], buf_values[k], first_row);
//        }
//        if((buf_j_idx[k] - first_row) >= mat_row){
//            printf("[%d] col Outof index for i=%d, j=%d, val=%lf, first row=%d\n", rank, buf_i_idx[k], buf_j_idx[k], buf_values[k], first_row);
//        }
//        printf("[%d] first row=%d, i=%d, x=%lf, val = %lf, result=%lf\n", rank, first_row, buf_i_idx[k], buf_x[buf_i_idx[k] - first_row], buf_values[k], buf_values[k] * buf_x[buf_i_idx[k] - first_row]);
//        printf("[%d] i=%d, j=%d, val=%lf, first row=%d\n", rank, buf_i_idx[k], buf_j_idx[k], buf_values[k], first_row);
        y[buf_i_idx[k] - first_row] += buf_values[k] * buf_x[buf_j_idx[k] - first_row];
    }
    return y;
}

int main(int argc, char *argv[]) {

    double t, comp_time, partition_time;
    int nprocs,     /* number of tasks/processes */
            rank;       /* id of task/process */

    int *buf_i_idx,     /* row index for all matrix elements */
            *buf_j_idx;     /* column index for all matrix elements */
    double *buf_values, /* value for all matrix elements */
            *buf_x;      /* value for all x vector elements */

    /*******************************************/

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int mat_size = 0, nonZero = 0, nonZeroPerRow = 0, total_run = 100, mat_row = 0, mat_col=0, first_row, last_row;

    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        mat_size = atoi(argv[1]);
        nonZeroPerRow = atoi(argv[2]);
        if (argc > 3)
            total_run = atoi(argv[3]);
    }

    mat_row = mat_size/nprocs;
    mat_col = mat_size;
    if(nonZeroPerRow > mat_row) {
        if(rank == MASTER) {
            printf("[%d] nonzero=%d, max nonzero=%d, number process=%d\n", rank, nonZeroPerRow, mat_row/2, nprocs);
        }
        nonZeroPerRow = mat_row /2;
    }
    nonZero = nonZeroPerRow * mat_row;
    if(nonZero <= 0){
        printf("[%d] Matrix can not be sized zero=%d\n", rank, nonZero);
    }
    first_row = rank * mat_row;
//    last_row = ((rank + 1) * mat_row) -1;
    buf_i_idx = (int *)malloc( nonZero * sizeof(int) );
    buf_j_idx = (int *)malloc( nonZero * sizeof(int) );
    buf_values = (double *)malloc( nonZero * sizeof(double));
    if(random_mat(buf_i_idx, buf_j_idx, buf_values, first_row, mat_row, nonZeroPerRow) != 1){
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nprocs, mat_size, nonZeroPerRow);
    }
    buf_x = (double *) malloc_or_exit(mat_row * sizeof(double));
    for (int i = 0; i < mat_row; i++) {
        buf_x[i] = 1.00;
    }
    /*printf("[%d] row buf = |", rank);
    for (int j = 0; j < nonZero; ++j) {
        printf("%d|",buf_i_idx[j]);
    }
    printf("\n");
    printf("[%d] col buf = |", rank);
    for (int j = 0; j < nonZero; ++j) {
        printf("%d|",buf_j_idx[j]);
    }
    printf("\n");
    printf("[%d] vlue buf = |", rank);
    for (int j = 0; j < nonZero; ++j) {
        printf("%lf|",buf_values[j]);
    }
    printf("\n");*/



    /* Matrix-vector multiplication for each processes */
    double totalTime = 0.0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    double *res;
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    for (int r = 0; r < total_run; ++r) {
        res = matMullComputationOnly(rank, buf_i_idx, buf_j_idx, buf_values, buf_x, first_row, mat_row, nonZero);
    }
    MPI_Barrier(MPI_COMM_WORLD);
//    printf("[%d] done multiplication\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
    totalTime = (MPI_Wtime() - t) * 1000.00;
    avg_time = totalTime / total_run;

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
        if ((checkFile = fopen("MPISpMVComputationResult.csv", "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen("MPISpMVComputationResult.csv", "a"))) {
                fprintf(stderr, "fopen: failed to open file MPISpMVComputationResult.csv");
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen("MPISpMVComputationResult.csv", "w"))) {
                fprintf(stderr, "fopen: failed to open file MPISpMVComputationResult.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "MatrixSize,MinTime,MaxTime,AvgTime,TotalRun,nProcess,NonZeroPerBlock,AvgOnDiagonalColumn\n");
        }

        fprintf(resultCSV, "%d,%10.3lf,%10.3lf,%10.3lf,%d,%d,%d,%d\n", mat_size, min_time, max_time, mean,
                total_run, nprocs, nonZero, nonZero);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file MPISpMVResult");
            exit(EXIT_FAILURE);
        }
    }

    free(buf_values);
    free(buf_i_idx);
    free(buf_j_idx);
    free(buf_x);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


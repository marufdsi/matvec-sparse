#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <jmorecfg.h>


#define TOTAL_RUNS 100

#define MAX_RANDOM_NUM (1<<20)
#define MASTER 0
#define EPSILON 1e-9

#include "mmio-wrapper.h"
#include "util.h"
#include "partition.h"

#include <omp.h>


double *matMull(int *row_ptr, int *col_ptr, double *val_ptr, double *x, int nRow, double *y) {

    /// multiplication
    #pragma omp parallel for schedule(dynamic, 64)
    for (int i = 0; i < nRow; ++i) {
        double tmp = 0;
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * x[col_ptr[k]];
//        y[i] = tmp;
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

    double comp_time = 0.0, avg_time = 0.0, *val_ptr, *x, *y;
    int total_run = 1000, skip=100, *row_ptr, *col_ptr, mat_row, _nnz, type=0;
    char *affinity, *inputFileName;
    int reserve_nodes =  omp_get_max_threads();
    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        inputFileName = argv[1];
        if (argc > 2)
            total_run = atoi(argv[2]);
        if (argc > 3)
            affinity = argv[3];
	if (argc > 4)
	  reserve_nodes = atoi(argv[4]);
    }

    if (read_coo_matrix_to_csr(inputFileName, &row_ptr, &col_ptr, &val_ptr, &mat_row, &_nnz) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }

    double nzPerRow = (double)_nnz/mat_row;

    y = (double *) calloc_or_exit(mat_row, sizeof(double));
    x = (double *) malloc_or_exit(mat_row * sizeof(double));
    for (int i = 0; i < mat_row; ++i) {
        x[i] = 1.0;
    }
    /// Start sparse matrix vector multiplication
    struct timespec start_time, end_time;
    for (int r = 0; r < total_run+skip; ++r) {
        if(r>=skip)
            clock_gettime(CLOCK_MONOTONIC, &start_time);
        // Multiplication
        matMull(row_ptr, col_ptr, val_ptr, x, mat_row, y);
        if(r>=skip) {
            clock_gettime(CLOCK_MONOTONIC, &end_time);
            comp_time += ((end_time.tv_sec * 1000 + (end_time.tv_nsec / 1.0e6)) -
                          (start_time.tv_sec * 1000 + (start_time.tv_nsec / 1.0e6)));
        }
    }
    avg_time = comp_time / total_run;

    int max_tid = omp_get_max_threads();
    char *ptr = strtok(inputFileName, "/");
    char *matrixName = strtok(strtok(NULL, "-"), ".");
    printf("Parallel SpMV computational time: %lf of matrix %s, size: %d, using %d threads\n", avg_time, matrixName, mat_row, max_tid);

    char outputFile[] = "OpenMP_Parallel_SpMV_on_Graphs.csv";
    FILE *resultCSV;
    FILE *checkFile;
    if ((checkFile = fopen(outputFile, "r")) != NULL) {
        // file exists
        fclose(checkFile);
        if (!(resultCSV = fopen(outputFile, "a"))) {
            fprintf(stderr, "fopen: failed to open %s file\n", outputFile);
            exit(EXIT_FAILURE);
        }
    } else {
        if (!(resultCSV = fopen(outputFile, "w"))) {
            fprintf(stderr, "fopen: failed to open file %s\n",outputFile);
            exit(EXIT_FAILURE);
        }
        fprintf(resultCSV, "GraphName,MatrixSize,AvgTime,TotalRun,Threads,NonZeroPerRow,NonZeroElements,AffinityType,ReserveNodes\n");
    }

    fprintf(resultCSV, "%s,%d,%10.3lf,%d,%d,%10.3lf,%d,%s,%d\n", matrixName, mat_row, avg_time, total_run, max_tid, nzPerRow, _nnz, affinity, reserve_nodes);
    if (fclose(resultCSV) != 0) {
        fprintf(stderr, "fopen: failed to open file %s\n", outputFile);
        exit(EXIT_FAILURE);
    }

    free(x);
    free(y);
    free(row_ptr);
    free(col_ptr);
    free(val_ptr);
    /* MPI: end */

    return 0;
}


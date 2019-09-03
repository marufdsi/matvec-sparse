#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <jmorecfg.h>

#include "mpi.h"

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
    #pragma omp for schedule(dynamic)
    for (int i = 0; i < nRow; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * x[col_ptr[k]];
    }
}

int createCSRMat(int **row_ptr, int **col_ptr, double **val_ptr, int mat_row, int _nnz){
    (*row_ptr) = (int *) calloc_or_exit(mat_row+1, sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(_nnz * sizeof(int));
    (*val_ptr) = (double *) malloc_or_exit(_nnz * sizeof(double));
    srand(time(NULL) * (omp_get_max_threads() + 1));
    int max_nnz_per_row = ceil((double)_nnz/mat_row);
    int **trackIndex = (int **) malloc_or_exit( mat_row * sizeof(int*));
    int *counter = (int *) calloc_or_exit(mat_row, sizeof(int));
    for (int i = 0; i < mat_row; ++i) {
        trackIndex[i] = (int *) malloc_or_exit(max_nnz_per_row * sizeof(int));
    }
    for (int i = 0; i < mat_row; ++i) {
        for (int j = 0; j < max_nnz_per_row; ++j) {
            trackIndex[i][j] = -1;
        }
    }
    int isExist = 0;
    for (int k = 0; k < _nnz; ++k) {
        int rowIdx = k%mat_row;
        int randColIdx;
        do {
            isExist = 0;
            randColIdx = rand() % mat_row;
            for (int i = 0; i < counter[rowIdx]; ++i) {
                if(trackIndex[rowIdx][i] == randColIdx) {
                    isExist = 1;
                    break;
                }
            }
        } while (isExist != 0);
        trackIndex[rowIdx][counter[rowIdx]++] = randColIdx;
        (*row_ptr)[rowIdx]++;
        (*col_ptr)[k] = randColIdx;
        (*val_ptr)[k] = 1.0;
    }
    for (int i = 0, cumsum = 0; i < mat_row; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = cumsum;
        cumsum += temp;
    }
    (*row_ptr)[mat_row] = _nnz;
    for (int i = 0; i < mat_row; ++i) {
        free(trackIndex[i]);
    }
    free(counter);
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

    double comp_time = 0.0, avg_time = 0.0, *val_ptr, *x, *y;
    int total_run = 1000, skip=100, *row_ptr, *col_ptr, mat_row, _nnz;


    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        mat_row = atoi(argv[1]);
        _nnz = atoi(argv[2]);
        if (argc > 3)
            total_run = atoi(argv[3]);
    }

    if (createCSRMat(&row_ptr, &col_ptr, &val_ptr, mat_row, _nnz) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }
    y = (double *) calloc_or_exit(mat_row, sizeof(double));
    x = (double *) malloc_or_exit(mat_row * sizeof(double));
    for (int i = 0; i < mat_row; ++i) {
        x[i] = 1.0;
    }
    /// Start sparse matrix vector multiplication
    struct timespec start_time, end_time;
    for (int r = 0; r < total_run+skip; ++r) {
        clock_gettime(CLOCK_MONOTONIC, &start_time);
        // Multiplication
        matMull(row_ptr, col_ptr, val_ptr, x, mat_row, y);
        clock_gettime(CLOCK_MONOTONIC, &end_time);
        comp_time += ((end_time.tv_sec * 1000 + (end_time.tv_nsec / 1.0e6)) - (start_time.tv_sec * 1000 + (start_time.tv_nsec / 1.0e6)));
    }
    avg_time = comp_time / total_run;

    int max_tid = omp_get_max_threads();
    printf("Parallel SpMV computational time: %lf of matrix size: %d, using %d threads\n", avg_time, mat_row, max_tid);

    FILE *resultCSV;
    FILE *checkFile;
    if ((checkFile = fopen("OMP_CSR_SpMV_Model.csv", "r")) != NULL) {
        // file exists
        fclose(checkFile);
        if (!(resultCSV = fopen("OMP_CSR_SpMV_Model.csv", "a"))) {
            fprintf(stderr, "fopen: failed to open file OMP_CSR_SpMV_Model.csv");
            exit(EXIT_FAILURE);
        }
    } else {
        if (!(resultCSV = fopen("OMP_CSR_SpMV_Model.csv", "w"))) {
            fprintf(stderr, "fopen: failed to open file OMP_CSR_SpMV_Model.csv");
            exit(EXIT_FAILURE);
        }
        fprintf(resultCSV, "MatrixSize,AvgTime,TotalRun,Threads,NonZeroPerRow,NonZeroPerBlock\n");
    }

    fprintf(resultCSV, "%d,%10.3lf,%d,%d,%10.3lf,%d\n", mat_row, avg_time, total_run, max_tid, (double)_nnz/mat_row, _nnz);
    if (fclose(resultCSV) != 0) {
        fprintf(stderr, "fopen: failed to open file OMP_CSR_SpMV_Model.csv");
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


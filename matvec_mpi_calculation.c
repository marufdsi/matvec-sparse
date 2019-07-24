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

/* Partition policy selection {EQUAL_NZ, EQUAL_ROWS} */
MPI_Datatype proc_info_type;
proc_info_t proc_info;
int first_row;

/*
 * Creates two MPI derived datatypes for the respective structs
 */
void create_mpi_datatypes(MPI_Datatype *proc_info_type) {
    MPI_Datatype oldtypes[2];
    MPI_Aint offsets[2];
    int blockcounts[2];

    /* create `proc_info_t` datatype */
    offsets[0] = 0;
    oldtypes[0] = MPI_INT;
    blockcounts[0] = 9;

    MPI_Type_create_struct(1, blockcounts, offsets, oldtypes, proc_info_type);
    MPI_Type_commit(proc_info_type);
}

double *matMullComputationOnly(int rank, int *buf_i_idx, int *buf_j_idx, double *buf_values, double *buf_x) {
    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(proc_info.M, sizeof(double));
    /// Sparse Matrix Vector Multiplication without Communication
    for (int k = 0; k < proc_info.NZ; k++) {
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

    create_mpi_datatypes(&proc_info_type);
    int mat_size = 0, nonZero = 0, total_run = 100;

    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        mat_size = atoi(argv[1]);
        nonZero = atoi(argv[2]);
        if (argc > 3)
            total_run = atoi(argv[3]);
    }

    proc_info.M = mat_size/nprocs;
    proc_info.N = mat_size;
    if(nonZero>(proc_info.M*proc_info.M))
        nonZero = proc_info.M * proc_info.M;
    proc_info.NZ = nonZero;
    proc_info.first_row = rank * proc_info.M;
    proc_info.last_row = (rank + 1) * proc_info.M;
    buf_i_idx = (int *)malloc( nonZero * sizeof(int) );
    buf_j_idx = (int *)malloc( nonZero * sizeof(int) );
    buf_values = (double *)malloc( nonZero * sizeof(double));
    if(random_mat(buf_i_idx, buf_j_idx, buf_values, proc_info.first_row, proc_info.last_row, proc_info.NZ) != 1)
    first_row = proc_info.first_row;
    buf_x = (double *) malloc_or_exit(proc_info.M * sizeof(double));
    for (int i = 0; i < proc_info.M; i++) {
        buf_x[i] = 1;
    }

    /* Matrix-vector multiplication for each processes */
    double totalTime = 0.0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    double *res;
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    for (int r = 0; r < total_run; ++r) {
        res = matMullComputationOnly(rank, buf_i_idx, buf_j_idx, buf_values, buf_x);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    totalTime = (MPI_Wtime() - t) * 1000.00;
    avg_time = totalTime / total_run;

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nprocs;

    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n",
               rank, min_time, max_time, mean, proc_info.NZ);
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

        fprintf(resultCSV, "%d,%10.3lf,%10.3lf,%10.3lf,%d,%d,%d,%d\n", proc_info.N, min_time, max_time, mean,
                total_run, nprocs, proc_info.NZ, proc_info.NZ);
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


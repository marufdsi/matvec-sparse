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
proc_info_t *proc_info;

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
    double *y = (double *) calloc_or_exit(proc_info[rank].M, sizeof(double));
    /// Sparse Matrix Vector Multiplication without Communication
    for (int k = 0; k < proc_info[rank].NZ; k++) {
        /*if (rank == 3){
            printf("[%d] i=%d, j=%d, val=%lf, first row=%d\n", rank, buf_i_idx[k], buf_j_idx[k], buf_values[k], proc_info[rank].first_row);
        }*/
        y[buf_i_idx[k] - proc_info[rank].first_row] += buf_values[k] * buf_x[buf_j_idx[k] - proc_info[rank].first_row];
    }
    return y;
}

int main(int argc, char *argv[]) {
    char *in_file, *out_file = NULL;

    double t, comp_time, partition_time;
    int nprocs,     /* number of tasks/processes */
            rank;       /* id of task/process */

    int *buf_i_idx,     /* row index for all matrix elements */
            *buf_j_idx;     /* column index for all matrix elements */
    double *buf_values, /* value for all matrix elements */
            *vec_x;      /* value for all x vector elements */

    /*******************************************/

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    create_mpi_datatypes(&proc_info_type);

    if (argc < 2 || argc > 3) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        in_file = argv[1];
        if (argc == 3)
            out_file = argv[2];
    }

    /// Initialize process info array
    proc_info = (proc_info_t *) malloc_or_exit(nprocs * sizeof(proc_info_t));
    /// Read input matrix
    if (rank_wise_read_matrix(in_file, &buf_i_idx, &buf_j_idx, &buf_values,
                              &proc_info[rank].M, &proc_info[rank].N, &proc_info[rank].NZ,
                              &proc_info[rank].first_row, &proc_info[rank].last_row, rank) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }
    vec_x = (double *) malloc_or_exit(proc_info[rank].N * sizeof(double));
    for (int i = 0; i < proc_info[rank].N; i++) {
        vec_x[i] = 1;
    }

    /// Share process info among all the processes
    MPI_Allgather(&proc_info[rank], 1, proc_info_type, proc_info, 1, proc_info_type, MPI_COMM_WORLD);

    /* Matrix-vector multiplication for each processes */
    double totalTime = 0.0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    double *res;
    for (int r = 0; r < TOTAL_RUNS; ++r) {
        MPI_Barrier(MPI_COMM_WORLD);
        t = MPI_Wtime();
        res = matMullComputationOnly(rank, buf_i_idx, buf_j_idx, buf_values, vec_x);
        totalTime += (MPI_Wtime() - t) * 1000.00;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    printf("[%d] Total run time: %lf\n", rank, totalTime);
    avg_time = totalTime / TOTAL_RUNS;
    /*printf("[%d] Results, y = |", rank);
    for (int i = 0; i<proc_info[rank].M; ++i) {
        printf("%lf|", res[i]);
    }
    printf("\n");*/
    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nprocs;

    int minNonZero = 0, maxNonZero = 0, avgNonZero = 0;
    MPI_Reduce(&proc_info[rank].NZ, &minNonZero, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&proc_info[rank].NZ, &maxNonZero, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&proc_info[rank].NZ, &avgNonZero, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    avgNonZero = avgNonZero / nprocs;

    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, Min NonZero: %d, Max NonZero: %d, Avg NonZero: %d\n",
               rank, min_time, max_time, mean, minNonZero, maxNonZero, avgNonZero);
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
                    "MatrixName,MinTime,MaxTime,AvgTime,TotalRun,nProcess,MinNonZero,MaxNonZero,AvgNonZero\n");
        }

        fprintf(resultCSV, "%s,%10.3lf,%10.3lf,%10.3lf,%d,%d,%d,%d,%d\n", in_file, min_time, max_time, mean, TOTAL_RUNS,
                nprocs, minNonZero, maxNonZero, avgNonZero);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file MPISpMVResult");
            exit(EXIT_FAILURE);
        }
    }

    free(buf_values);
    free(buf_i_idx);
    free(buf_j_idx);
    free(vec_x);
    MPI_Finalize();

    return 0;
}


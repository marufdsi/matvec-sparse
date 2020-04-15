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

MPI_Datatype procs_info_type;
enum tag {
    REQUEST_TAG, RECEIVE_TAG
};

//typedef float f_type;

f_type *matMull(int rank, int *row_ptr, int *col_ptr, f_type *val_ptr, f_type *x, int nRow, int startCol, f_type *y) {

    /// multiplication
    for (int i = 0; i < nRow; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * x[col_ptr[k] - startCol];
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
    double comp_time = 0, bcast_time = 0.0, matmul_time = 0.0, reduce_time = 0.0, min_time = 0.0, max_time = 0.0,
            avg_time = 0.0, mean = 0.0, avg_bcast_time = 0.0, avg_matmul_time = 0.0, avg_reduce_time = 0.0;
    int total_run = 30, skip = 5, nRanks, rank, knl = 0, TOTAL_MAT_MUL = 20, nodes = 0;

    int *row_ptr, *col_ptr;
    f_type *val_ptr, *x, *y;
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
        in_file = argv[1];
        if (argc > 2)
            nodes = atoi(argv[2]);
        if (argc > 3)
            knl = atoi(argv[3]);
    }
    if (csr_read_2D_partitioned_mat(in_file, &row_ptr, &col_ptr, &val_ptr, &ranks_info, sqrRank, rank) != 0) {
        printf("Error in the process=%d\n", rank);
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }
    printf("[%d] Done Reading", rank);
    printf(" M(%d)=%d !\n", rank, ranks_info[rank].M);
    y = (f_type *) malloc_or_exit(ranks_info[rank].M * sizeof(f_type));
    x = (f_type *) malloc_or_exit(ranks_info[rank].M * sizeof(f_type));
    for (int i = 0; i < ranks_info[rank].M; ++i) {
        x[i] = 1.0;
        y[i] = 0.0;
    }
    printf("[%d] Done Initialization!\n", rank);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank == MASTER){
        printf("Ready to perform mat-mul\n");
    }
    MPI_Finalize();

    return 0;
//    printf("[%d] Vector creation done!\n", rank);

    /// Share process info among all the processes
    MPI_Allgather(&ranks_info[rank], 1, procs_info_type, procs_info, 1, procs_info_type, MPI_COMM_WORLD);
    /// Start sparse matrix vector multiplication for each rank
    MPI_Barrier(MPI_COMM_WORLD);
    struct timespec start, end, b_start, b_end, r_start, r_end, m_start, m_end;
    for (int r = 0; r < total_run + skip; ++r) {
        for (int mul = 0; mul < TOTAL_MAT_MUL; ++mul) {
            clock_gettime(CLOCK_MONOTONIC, &start);
            clock_gettime(CLOCK_MONOTONIC, &b_start);
            //broadcast X along column communicator
            MPI_Bcast(x, ranks_info[rank].M, MPI_FLOAT, col_rank,
                      commcol); //col_rank is the one with the correct information
            if (r >= skip) {
                clock_gettime(CLOCK_MONOTONIC, &b_end);
                bcast_time += ((b_end.tv_sec * 1000 + (b_end.tv_nsec / 1.0e6)) -
                               (b_start.tv_sec * 1000 + (b_start.tv_nsec / 1.0e6)));
            }


            clock_gettime(CLOCK_MONOTONIC, &m_start);
            // Multiplication
            matMull(rank, row_ptr, col_ptr, val_ptr, x, ranks_info[rank].M, col_rank * ranks_info[rank].M, y);
            if (r >= skip) {
                clock_gettime(CLOCK_MONOTONIC, &m_end);
                matmul_time += ((m_end.tv_sec * 1000 + (m_end.tv_nsec / 1.0e6)) -
                                (m_start.tv_sec * 1000 + (m_start.tv_nsec / 1.0e6)));
            }

            clock_gettime(CLOCK_MONOTONIC, &r_start);
            //reduce Y along row communicator
            MPI_Reduce(y, x, ranks_info[rank].M, MPI_FLOAT, MPI_SUM, row_rank, commrow);
            if (r >= skip) {
                clock_gettime(CLOCK_MONOTONIC, &r_end);
                clock_gettime(CLOCK_MONOTONIC, &end);
                reduce_time += ((r_end.tv_sec * 1000 + (r_end.tv_nsec / 1.0e6)) -
                                (r_start.tv_sec * 1000 + (r_start.tv_nsec / 1.0e6)));
                comp_time += ((end.tv_sec * 1000 + (end.tv_nsec / 1.0e6)) -
                              (start.tv_sec * 1000 + (start.tv_nsec / 1.0e6)));
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
//    comp_time = (MPI_Wtime() - start_time) * 1000.00;
    avg_time = comp_time / total_run;
    avg_bcast_time = bcast_time / total_run;
    avg_matmul_time = matmul_time / total_run;
    avg_reduce_time = reduce_time / total_run;

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nRanks;
    MPI_Comm_free(&commrow);
    MPI_Comm_free(&commcol);
    int max_nnz = 0, sum_nnz = 0, avg_row = 0;
    double max_nnz_per_row = 0.0, avg_nnz_per_row = 0.0, nnz_per_row = (double) procs_info[rank].NZ / procs_info[rank].M,
    avg_nnz = 0.0, sd = 0.0, avg_sd = 0.0;

    for (int j = 1; j < procs_info[rank].M+1; ++j) {
        sd += (row_ptr[j] - row_ptr[j-1] - nnz_per_row) * (row_ptr[j] - row_ptr[j-1] - nnz_per_row);
    }
    sd = sd/procs_info[rank].M;
    sd = sqrt(sd);
    MPI_Reduce(&procs_info[rank].NZ, &max_nnz, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&procs_info[rank].NZ, &sum_nnz, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&procs_info[rank].M, &avg_row, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&nnz_per_row, &avg_nnz_per_row, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&nnz_per_row, &max_nnz_per_row, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&sd, &avg_sd, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
    avg_nnz = sum_nnz / nRanks;
    avg_nnz_per_row = avg_nnz_per_row/nRanks;
    avg_sd = avg_sd/nRanks;
    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n",
               rank, min_time, max_time, mean, procs_info[rank].NZ);
//        char *_ptr = strtok(in_file, "/");
        char * ptr = strtok(in_file, "/");
        char *file[1025];
        int f_i=0;
        while(ptr != NULL)
        {
            file[f_i++] = ptr;
            ptr = strtok(NULL, "/");
        }
        ptr = strtok(file[f_i-1], ".");
        char *matrixName = ptr;
        char outputFile[100] = "Skylake_CSR_Random_BrCast_Reduce_SpMV.csv";
        if (knl > 0)
            strcpy(outputFile, "KNL_CSR_Random_BrCast_Reduce_SpMV.csv");

        FILE *resultCSV;
        FILE *checkFile;
        if ((checkFile = fopen(outputFile, "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen(outputFile, "a"))) {
                fprintf(stderr, "fopen: failed to open file %s", outputFile);
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen(outputFile, "w"))) {
                fprintf(stderr, "fopen: failed to open file %s", outputFile);
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "Name,MatrixSize,AvgRow,MinTime,MaxTime,AvgTime,AvgBcastTime,AvgMatmulTime,AvgReduceTime,TotalRun,nProcess,NonZeroPerRow,MaxNonZeroPerRow,AvgNonZeroPerBlock,MaxNonZeroPerBlock,Nodes,AvgNPRSD,DataType\n");
        }

        fprintf(resultCSV, "%s,%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%d,%d,%lf,%lf,%lf,%d,%d,%lf,%d\n", matrixName, procs_info[rank].N,
                avg_row, min_time, max_time,
                mean, avg_bcast_time, avg_matmul_time, avg_reduce_time, total_run, nRanks, avg_nnz_per_row, max_nnz_per_row, avg_nnz,
                max_nnz, nodes, avg_sd, sizeof(f_type));
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file %s", outputFile);
            exit(EXIT_FAILURE);
        }
    }
    free(x);
    free(y);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


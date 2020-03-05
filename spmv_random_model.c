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
#ifndef DIAGONAL_MATRIX
#define DIAGONAL_MATRIX true
#endif

#ifndef RANDOM_MATRIX
#define RANDOM_MATRIX false
#endif

#ifndef CSR_MATRIX
#define CSR_MATRIX false
#endif

double *matMull(int rank, int *row_ptr, int *col_ptr, double *val_ptr, double *x, int nRow, int startCol, double *y) {

    /// multiplication
    for (int i = 0; i < nRow; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i+1]; ++k) {
//            if((col_ptr[k] - startCol) < 0)
//                printf("Negative index found\n");
//            if((col_ptr[k] - startCol) >= nRow)
//                printf("index out of bound\n");
            y[i] += val_ptr[k] * x[col_ptr[k] - startCol];
        }
    }
}

int  create_random_matrix(int **row_ptr, int **col_ptr, double **val_ptr, int m, int n, int nnz_per_row, int startCol, int rank){
    (*row_ptr) = (int *) calloc_or_exit(m+1, sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(m*nnz_per_row * sizeof(int));
    (*val_ptr) = (double *) malloc_or_exit(m*nnz_per_row * sizeof(double));
    srand(time(NULL) * (rank + 1));

    int *trackIndex, idx = 0;
    for (int i = 0; i < m; ++i) {
        trackIndex = (int *)calloc_or_exit( n, sizeof(int));
        for (int j = 0; j < nnz_per_row; ++j) {
            int randColIdx;
            do {
                randColIdx = rand() % n;
            } while (trackIndex[randColIdx] != 0);
            trackIndex[ randColIdx] = 1.0;
            (*row_ptr)[i]++;
            (*col_ptr)[idx] = startCol + randColIdx;
            (*val_ptr)[idx] = ((double)(randColIdx%10) +1);
            idx++;
        }
        free(trackIndex);
    }
    return 0;
}

int  create_random_diagonal_matrix(int **row_ptr, int **col_ptr, double **val_ptr, int m, int nnz_per_row, int startCol, int rank){
    (*row_ptr) = (int *) calloc_or_exit(m+1, sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(m*nnz_per_row * sizeof(int));
    (*val_ptr) = (double *) malloc_or_exit(m*nnz_per_row * sizeof(double));
    srand(time(NULL) * (rank + 1));

    int *trackIndex, idx = 0;
    trackIndex = (int *)malloc_or_exit( m * sizeof(int));
    for (int i = 0; i < m; ++i) {
        for(int l=0; l<m; ++l)
            trackIndex[l] = 0;
        for (int j = 0; j < nnz_per_row; ++j) {
            int randColIdx;
            do {
                randColIdx = rand() % m;
            } while (trackIndex[randColIdx] != 0);
            trackIndex[ randColIdx] = 1.0;
            (*row_ptr)[i]++;
            (*col_ptr)[idx] = startCol + randColIdx;
            (*val_ptr)[idx] = ((double)(randColIdx%10) +1);
            idx++;
        }
    }
    for(int i=0, prev=0; i<m+1; ++i){
        int new = (*row_ptr)[i] + prev;
        (*row_ptr)[i] = prev;
        prev = new;
    }
    free(trackIndex);
    return 0;
}

int createCSRMat(int **row_ptr, int **col_ptr, double **val_ptr, int mat_row, int nnz_per_block, int startCol, int rank){
    (*row_ptr) = (int *) calloc_or_exit(mat_row+1, sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(nnz_per_block * sizeof(int));
    (*val_ptr) = (double *) malloc_or_exit(nnz_per_block * sizeof(double));
    srand(time(NULL) * (rank + 1));
    int max_nnz_per_row = ceil((double)nnz_per_block/mat_row);
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
    for (int k = 0; k < nnz_per_block; ++k) {
        int rowIdx = rank%2 == 0 ? (mat_row-1)-(k%mat_row) : (k%mat_row);
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
        (*col_ptr)[k] = startCol + randColIdx;
        (*val_ptr)[k] = 1.0;
    }
    for (int i = 0, cumsum = 0; i < mat_row; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = cumsum;
        cumsum += temp;
    }
    (*row_ptr)[mat_row] = nnz_per_block;
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

    char *in_file;
    double comp_time = 0.0, bcast_time = 0.0, matmul_time = 0.0, reduce_time = 0.0, min_time = 0.0, max_time = 0.0,
            avg_time = 0.0, mean = 0.0, avg_bcast_time = 0.0, avg_matmul_time = 0.0, avg_reduce_time = 0.0,
            *val_ptr, *x, *y, *y_seq;
    int total_run = 30, skip=5, nRanks, rank, *row_ptr, *col_ptr, _size, mat_row, nnz_per_block, TOTAL_MAT_MUL = 20,
            nodes = 0, procs_per_node = 0;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        _size = atoi(argv[1]);
        mat_row = atoi(argv[2]);
        nnz_per_block = atoi(argv[3]);
        if (argc > 4)
            total_run = atoi(argv[4]);
        if (argc > 5)
            nodes = atoi(argv[5]);
        if (argc > 6)
            procs_per_node = atoi(argv[6]);
    }

    if(rank == MASTER){
        printf("matrix ize: %d, row: %d, nzpb: %d, run: %d, nodes: %d, ppn: %d\n", _size, mat_row, nnz_per_block,
                total_run, nodes, procs_per_node);
    }
    int sqrRank = sqrt(nRanks);
    int row_rank = rank / sqrRank; //which col of proc am I
    int col_rank = rank % sqrRank; //which row of proc am I

    //initialize communicators
    MPI_Comm commrow;
    MPI_Comm_split(MPI_COMM_WORLD, row_rank, rank, &commrow);

    MPI_Comm commcol;
    MPI_Comm_split(MPI_COMM_WORLD, col_rank, rank, &commcol);

#if CSR_MATRIX
    printf("CSR matrix called\n");
    if (createCSRMat(&row_ptr, &col_ptr, &val_ptr, mat_row, nnz_per_block, col_rank * mat_row, rank) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }
#endif
//#if DIAGONAL_MATRIX
    printf("Diagonal matrix called\n");
    if (create_random_diagonal_matrix(&row_ptr, &col_ptr, &val_ptr, mat_row, nnz_per_block/mat_row, col_rank * mat_row, rank) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }
//#endif
    y = (double *) malloc_or_exit(mat_row * sizeof(double));
    x = (double *) malloc_or_exit(mat_row * sizeof(double));
    for (int i = 0; i < mat_row; ++i) {
        x[i] = 1.0;
    }
    if (rank == MASTER){
        printf("[%d] Matrix creation done\n", rank);
    }
    /// Start sparse matrix vector multiplication for each rank
    MPI_Barrier(MPI_COMM_WORLD);
    struct timespec start, end, b_start, b_end, r_start, r_end, m_start, m_end;
    for (int r = 0; r < total_run+skip; ++r) {
        for (int mul = 0; mul < TOTAL_MAT_MUL; ++mul) {
            for (int i = 0; i < mat_row; ++i) {
                y[i] = 0.0;
            }
            clock_gettime(CLOCK_MONOTONIC, &start);
            clock_gettime(CLOCK_MONOTONIC, &b_start);
            //broadcast X along column communicator
            MPI_Bcast(x, mat_row, MPI_FLOAT, col_rank, commcol); //col_rank is the one with the correct information
            if (r >= skip) {
                clock_gettime(CLOCK_MONOTONIC, &b_end);
                bcast_time += ((b_end.tv_sec * 1000 + (b_end.tv_nsec / 1.0e6)) -
                        (b_start.tv_sec * 1000 + (b_start.tv_nsec / 1.0e6)));
            }

            clock_gettime(CLOCK_MONOTONIC, &m_start);
            // Multiplication
            matMull(rank, row_ptr, col_ptr, val_ptr, x, mat_row, col_rank * mat_row, y);
            if (r >= skip) {
                clock_gettime(CLOCK_MONOTONIC, &m_end);
                matmul_time += ((m_end.tv_sec * 1000 + (m_end.tv_nsec / 1.0e6)) -
                                (m_start.tv_sec * 1000 + (m_start.tv_nsec / 1.0e6)));
            }

            clock_gettime(CLOCK_MONOTONIC, &r_start);
            //reduce Y along row communicator
            MPI_Reduce(y, x, mat_row, MPI_FLOAT, MPI_SUM, row_rank, commrow);
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
//        char outputFile[100] = "CSR_SpMV_Model_of_Random_BrCast_Reduce.csv";
        char outputFile[100] = "Skylake_CSR_SpMV_Model_of_Diagonal_BrCast_Reduce.csv";

//        char outputFile[100] = "Skylake_CSR_SpMV_Model_of_Random_BrCast_Reduce.csv";

//        char outputFile[100] = "CSR_SpMV_Model_of_Random_BrCast_Reduce.csv";
        FILE *resultCSV;
        FILE *checkFile;
        if ((checkFile = fopen(outputFile, "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen(outputFile, "a"))) {
                fprintf(stderr, "fopen: failed to open file %s\n", outputFile);
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen(outputFile, "w"))) {
                fprintf(stderr, "fopen: failed to open file %s\n", outputFile);
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "MatrixSize,PartitionRow,MinTime,MaxTime,AvgTime,AvgBcastTime,AvgMatmulTime,AvgReduceTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock,Nodes,ProcessPerNode\n");
        }

        fprintf(resultCSV, "%d,%d,%10.3lf,%10.3lf,%10.3lf,%10.3lf,%10.3lf,%10.3lf,%d,%d,%10.3lf,%d,%d,%d\n", _size, mat_row,
                min_time, max_time, mean, avg_bcast_time, avg_matmul_time, avg_reduce_time, total_run, nRanks,
                nnz_per_row, nnz_per_block, nodes, procs_per_node);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file %s\n", outputFile);
            exit(EXIT_FAILURE);
        }
    }
    free(x);
    free(y);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


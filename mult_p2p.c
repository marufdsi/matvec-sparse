#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include<time.h>
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

void localMatMull(int *row_ptr, int *col_ptr, double *val_ptr, double *buf_x, int first_row, int mat_row, double *y) {
    /// Local elements multiplication
    for (int i = 0; i < mat_row; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * buf_x[col_ptr[k] - first_row];
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

    char *in_file, *out_file = NULL;
    double comp_time = 0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    int total_run = 100, skip=100, nRanks, rank, knl=0;
    int *row_ptr, *on_diagonal_row, *col_ptr, *on_diagonal_col;
    double *val_ptr, *on_diagonal_val, *buf_x, *res;
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

    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        in_file = argv[1];
        if (argc > 2)
            knl = atoi(argv[2]);
    }
    int offDiagonalElements = 0;
    if (rank_wise_read_matrix_csr(in_file, &row_ptr, &col_ptr, &val_ptr, &ranks_info, rank,
                                  &offDiagonalElements) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }

    if (ranks_info[rank].NZ <= 0) {
        printf("[%d] Matrix can not be sized zero=%d\n", rank, ranks_info[rank].NZ);
        return 0;
    }

    int diagonal_elements = ranks_info[rank].NZ - offDiagonalElements;
    if(diagonal_elements>0) {
        on_diagonal_row = (int *) malloc_or_exit((ranks_info[rank].M + 1) * sizeof(int));
        on_diagonal_col = (int *) malloc_or_exit(diagonal_elements * sizeof(int));
        on_diagonal_val = (double *) malloc_or_exit(diagonal_elements * sizeof(double));
        on_diagonal_row[0] = 0;
    }
    int on_diag_idx = 0;
    int in_diagonal_bandwidth = 0, bandwidth = 0;
    for (int k = 0; k < ranks_info[rank].M; ++k) {
        int l_col =ranks_info[rank].last_row, h_col=ranks_info[rank].first_row, max_col =ranks_info[rank].last_row, min_col=ranks_info[rank].first_row;
        for (int l = row_ptr[k]; l < row_ptr[k + 1]; ++l) {
            if (in_diagonal(col_ptr[l], ranks_info[rank].first_row, ranks_info[rank].last_row)) {
                on_diagonal_col[on_diag_idx] = col_ptr[l];
                on_diagonal_val[on_diag_idx] = val_ptr[l];
                on_diag_idx++;
                if(col_ptr[l]<l_col){
                    l_col = col_ptr[l];
                }
                if(col_ptr[l]>h_col){
                    h_col = col_ptr[l];
                }
            }
            if(col_ptr[l]<min_col){
                min_col = col_ptr[l];
            }
            if(col_ptr[l]>max_col){
                max_col = col_ptr[l];
            }
        }
        if(in_diagonal_bandwidth < (h_col - l_col + 1)){
            in_diagonal_bandwidth = (h_col - l_col + 1);
        }
        if(bandwidth < (max_col - min_col + 1)){
            bandwidth = (max_col - min_col + 1);
        }
        if(diagonal_elements>0)
            on_diagonal_row[k + 1] = on_diag_idx;
    }
    free(row_ptr);
    free(col_ptr);
    free(val_ptr);

    /// Share process info among all the processes
    MPI_Allgather(&ranks_info[rank], 1, procs_info_type, procs_info, 1, procs_info_type, MPI_COMM_WORLD);
    int first_row = procs_info[rank].first_row, mat_row = procs_info[rank].M;
    /// Create vector x and fill with 1.0
    buf_x = (double *) malloc_or_exit(mat_row* sizeof(double));
    for (int i = 0; i < mat_row; i++) {
        buf_x[i] = 1.00;
    }

    /// allocate memory for vectors and submatrixes
    double *y = (double *) calloc_or_exit(mat_row, sizeof(double));
    /// Start sparse matrix vector multiplication for each rank
    struct timespec start, end;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int r = 0; r < total_run+skip; ++r) {
        clock_gettime(CLOCK_MONOTONIC, &start);
        localMatMull(on_diagonal_row, on_diagonal_col, on_diagonal_val, buf_x, first_row, mat_row, y);
        if(r>=skip) {
            clock_gettime(CLOCK_MONOTONIC, &end);
            comp_time += ((end.tv_sec * 1000 + (end.tv_nsec / 1.0e6)) - (start.tv_sec * 1000 + (start.tv_nsec / 1.0e6)));
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    avg_time = comp_time / total_run;

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nRanks;

    int max_nnz = 0, max_row = 0, max_in_diag_band_width = 0, max_band_width = 0;
    double max_nnz_per_row = 0.0, nnz_per_row = procs_info[rank].NZ/procs_info[rank].M;
    MPI_Reduce(&procs_info[rank].NZ, &max_nnz, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&procs_info[rank].M, &max_row, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&in_diagonal_bandwidth, &max_in_diag_band_width, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&bandwidth, &max_band_width, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&nnz_per_row, &max_nnz_per_row, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);

//    char *outputFIle = (char *) malloc_or_exit(100 * sizeof(char));
//    strcpy(outputFIle, "CSR_SpMV_Model.csv");
//    strcpy(outputFIle, "CSR_SpMV_Model_Diagonal.csv");
//    strcpy(outputFIle, "CSR_SpMV_Model_Random.csv");
    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n",
               rank, min_time, max_time, mean, procs_info[rank].NZ);
        char *_ptr = strtok(in_file, "/");
        char *matrixName = strtok(strtok(NULL, "-"), ".");
        char outputFile[100] = "Skylake_CSR_Local_Matmul_on_MPI.csv";
        if(knl > 0)
            strcpy(outputFile, "KNL_CSR_Local_Matmul_on_MPI.csv");

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
                    "Name,MatrixSize,MaxRow,MinTime,MaxTime,AvgTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock,DiagonalElem,OffDiagonalElem,Bandwidth,InDiagonalBandwidth\n");
        }

        fprintf(resultCSV, "%s,%d,%d,%lf,%lf,%lf,%d,%d,%lf,%d,%d,%d,%d,%d\n", matrixName,procs_info[rank].N, max_row, min_time, max_time,
                mean, total_run, nRanks, max_nnz_per_row, max_nnz, diagonal_elements, offDiagonalElements, max_band_width, max_in_diag_band_width);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file %s", outputFile);
            exit(EXIT_FAILURE);
        }
    }

    if(diagonal_elements>0) {
        free(on_diagonal_row);
        free(on_diagonal_col);
        free(on_diagonal_val);
    }
    free(buf_x);
    free(y);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


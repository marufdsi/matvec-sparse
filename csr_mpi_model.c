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

MPI_Datatype proc_info_type;
enum tag {
    REQUEST_TAG, RECEIVE_TAG
};

struct Key {
    int col;
};
struct Value {
    double val;
};
struct Map {
    struct Key key;
    struct Value value;
};

double getVal(struct Map *map, int col) {
    while (map != NULL) {
        if ((*map).key.col == col) {
            return (*map).value.val;
        }
        map++;
    }
    printf("Error!!! column=%d not found\n", col);
    return 0;
}

void testMap(int rank) {
    if (rank == MASTER) {
        struct Map *map;
        map = (struct Map *) malloc_or_exit(10 * sizeof(struct Map));
        for (int i = 0; i < 10; ++i) {
            map[i].key.col = i * 10;
            map[i].value.val = (double) (i * 500);
        }

        for (int k = 0; k < 10; ++k) {
            printf("[%d] Map key=%d, value=%lf\n", k, map[k].key.col, map[k].value.val);
        }
        printf("Get value of col=%d from map=%lf\n", 10, getVal(map, 10));
    }
}

int getRank(int nRanks, proc_info_t *ranks_info, int column){
    for (int r = 0; r < nRanks; ++r) {
        if (column >= ranks_info[r].first_row && column <= ranks_info[r].last_row)
            return r;
    }
    printf("Error!! %d Column does not belong to any ranks\n", column);
    return -1;
}

double *matMull(int rank, proc_info_t *ranks_info, int nRanks, int *row_ptr, int *col_ptr, double *val_ptr, double *buf_x,
        int mat_row, int **send_col_idx, int *perRankDataRecv, int **reqColFromRank, int *perRankDataSend, int nColRecv) {

    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(mat_row, sizeof(double));
/// receiving blocks storage
    double **recv_buf = (double **) malloc_or_exit(nRanks * sizeof(double));
    double **recvColFromRanks = (int *) malloc_or_exit(nRanks * sizeof(double));
    for (int r = 0; r < nRanks; ++r){
        if (perRankDataRecv[r] > 0) {
            recv_buf[r] = (double *) malloc_or_exit(perRankDataRecv[r] * sizeof(double));
            recvColFromRanks[r] = (double *) malloc_or_exit((ranks_info[r].last_row - ranks_info[r].first_row + 1) * sizeof(double));
        }
    }
    /// MPI request storage
    MPI_Request *send_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
    MPI_Request *recv_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
    int reqMade = 0;
    for (int r = 0; r < nRanks; ++r) {
        if (r == rank || perRankDataRecv[r] == 0) {
            recv_reqs[r] = MPI_REQUEST_NULL;
            continue;
        }
        reqMade++;
        /// Receive the block (when it comes)
        MPI_Irecv(recv_buf[r], perRankDataRecv[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &recv_reqs[r]);
    }

    /// Reply to the requests.
    double **send_buf_data = (double **) malloc_or_exit(nRanks * sizeof(double));
    for (int r = 0; r < nRanks; ++r) {
        if (perRankDataSend[r] > 0) {
            send_buf_data[r] = (double *) malloc_or_exit(perRankDataSend[r] * sizeof(double));
            for (int i = 0; i < perRankDataSend[r]; ++i)
                send_buf_data[r][i] = buf_x[send_col_idx[r][i]];
            MPI_Isend(send_buf_data[r], perRankDataSend[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &send_reqs[r]);
        } else
            send_reqs[r] = MPI_REQUEST_NULL;
    }

    /// Local elements multiplication
    for (int i = 0; i < mat_row; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            if (in_diagonal(col_ptr[k], ranks_info[rank].first_row, ranks_info[rank].last_row))
                y[i] += val_ptr[k] * buf_x[col_ptr[k]];
        }
    }

    /// need to update the initialization
//    struct Map *map = (struct Map *) malloc_or_exit(nColRecv * sizeof(struct Map));
    int r, index = 0;
    for (int q = 0; q < reqMade; q++) {
        MPI_Waitany(nRanks, recv_reqs, &r, MPI_STATUS_IGNORE);
        assert(r != MPI_UNDEFINED);

        /// fill x array with new elements.
        for (int i = 0; i < perRankDataRecv[r]; i++) {
//            map[index].key.col = reqColFromRank[r][i];
//            map[index].value.val = recv_buf[r][i];
//            index++;
            recvColFromRanks[r][reqColFromRank[r][i] - ranks_info[r].first_row] = recv_buf[r][i];
        }
    }

    /// Global elements multiplication
    for (int i = 0; i < mat_row; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            if (!in_diagonal(col_ptr[k], ranks_info[rank].first_row, ranks_info[rank].last_row)) {
                int r = getRank(nRanks, ranks_info, col_ptr[k]);
                y[i] += val_ptr[k] * recvColFromRanks[r][col_ptr[k] - ranks_info[r].first_row];
            }
        }
    }

    /// Wait until send request delivered to through network.
    MPI_Waitall(nRanks, send_reqs, MPI_STATUS_IGNORE);
    for (int r = 0; r < nRanks; r++) {
        if (perRankDataRecv[r] > 0) {
            free(recv_buf[r]);
            free(recvColFromRanks[r]);
        }
        if (perRankDataSend[r] > 0)
            free(send_buf_data[r]);
    }

    free(recv_buf);
    free(send_buf_data);
    free(recvColFromRanks);

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

    proc_info_t *proc_info;
    /*******************************************/

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    testMap(rank);
    MPI_Finalize();
    return 0;

    int mat_size = 0, nonZero = 0, nonZeroPerRow = 0, total_run = 100, mat_row = 0, mat_col = 0;

    if (argc < 3) {
        printf("Usage: %s matrix_row matrix_col nonZeroPerRow [Total_Runs]\n", argv[0]);
        return 0;
    } else {
        mat_row = atoi(argv[1]);
        mat_col = atoi(argv[2]);
        nonZeroPerRow = atoi(argv[3]);
        if (argc > 3)
            total_run = atoi(argv[4]);
    }

    mat_size = mat_row * nprocs;
    nonZero = nonZeroPerRow * mat_row;
    if (nonZeroPerRow <= 0) {
        printf("[%d], There will must one non zero column in the matrix in every row\n");
        return 0;
    }
    if (nonZeroPerRow > mat_row) {
        if (rank == MASTER) {
            printf("[%d] nonzero=%d, max nonzero=%d, number process=%d\n", rank, nonZeroPerRow, mat_row, nprocs);
        }
        nonZeroPerRow = mat_row;
        nonZero = nonZeroPerRow * mat_row;
    }
    if (nonZero <= 0) {
        printf("[%d] Matrix can not be sized zero=%d\n", rank, nonZero);
    }

    row_ptr = (int *) malloc((mat_row + 1) * sizeof(int));
    col_ptr = (int *) malloc(nonZero * sizeof(int));
    val_ptr = (double *) malloc(nonZero * sizeof(double));
    if (csr_random_mat(row_ptr, col_ptr, val_ptr, mat_row, nonZeroPerRow) != 1) {
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nprocs, mat_size,
               nonZeroPerRow);
    }
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
        res = matMull(rank, proc_info, row_ptr, col_ptr, val_ptr, buf_x, mat_row);
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
               rank, min_time, max_time, mean, nonZero);
        FILE *resultCSV;
        FILE *checkFile;
        if ((checkFile = fopen("CSR_MPI_SpMV.csv", "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen("CSR_MPI_SpMV.csv", "a"))) {
                fprintf(stderr, "fopen: failed to open file CSR_MPI_SpMV.csv");
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen("CSR_MPI_SpMV.csv", "w"))) {
                fprintf(stderr, "fopen: failed to open file CSR_MPI_SpMV.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "MatrixSize,MinTime,MaxTime,AvgTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock\n");
        }

        fprintf(resultCSV, "%d,%10.3lf,%10.3lf,%10.3lf,%d,%d,%d,%d\n", mat_size, min_time, max_time, mean,
                total_run, nprocs, nonZeroPerRow, nonZero);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file MPISpMVResult");
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


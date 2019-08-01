#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

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

void testMap(int rank) {
    if (rank == MASTER) {
        Map *map;
        map = (Map *) malloc_or_exit(10 * sizeof(Map));
        for (int i = 0; i < 10; ++i) {
            map[i].key.col = i * 10;
            map[i].value.val = (double) (i * 500);
        }

        for (int k = 0; k < 10; ++k) {
            printf("[%d] Map key=%d, value=%lf\n", k, map[k].key.col, map[k].value.val);
        }
        printf("Get value of col=%d from map=%lf\n", 10, getVal(map, 10, 0));
    }
}

int getRank(int nRanks, proc_info_t *procs_info, int column) {
    for (int r = 0; r < nRanks; ++r) {
        if (column >= procs_info[r].first_row && column <= procs_info[r].last_row)
            return r;
    }
    printf("Error!! %d Column does not belong to any ranks\n", column);
    return -1;
}

double *
matMull(int rank, proc_info_t *procs_info, int nRanks, int *row_ptr, int *col_ptr, double *val_ptr, double *buf_x,
        int **send_col_idx, int *perRankDataRecv, int **reqColFromRank, int *perRankDataSend, int reqRequired,
        int nRanksExpectCol) {

    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(procs_info[rank].M, sizeof(double));
/// receiving blocks storage
    double **recv_buf, **recvColFromRanks;
    MPI_Request *send_reqs, *recv_reqs;
    int reqMade = 0;
    if (reqRequired > 0) {
        recv_buf = (double **) malloc_or_exit(nRanks * sizeof(double));
        recvColFromRanks = (double **) malloc_or_exit(nRanks * sizeof(double));
        for (int r = 0; r < nRanks; ++r) {
            if (perRankDataRecv[r] > 0) {
                recv_buf[r] = (double *) malloc_or_exit(perRankDataRecv[r] * sizeof(double));
                recvColFromRanks[r] = (double *) malloc_or_exit(
                        (procs_info[r].last_row - procs_info[r].first_row + 1) * sizeof(double));
            }
        }
        /// MPI request storage
        recv_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
        for (int r = 0; r < nRanks; ++r) {
            if (r == rank || perRankDataRecv[r] == 0) {
                recv_reqs[r] = MPI_REQUEST_NULL;
                continue;
            }
            reqMade++;
            /// Receive the block (when it comes)
            MPI_Irecv(recv_buf[r], perRankDataRecv[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &recv_reqs[r]);
        }
    }

    double **send_buf_data;
    if (nRanksExpectCol > 0) {
        send_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
        /// Reply to the requests.
        send_buf_data = (double **) malloc_or_exit(nRanks * sizeof(double));
        for (int r = 0; r < nRanks; ++r) {
            if (perRankDataSend[r] > 0) {
                send_buf_data[r] = (double *) malloc_or_exit(perRankDataSend[r] * sizeof(double));
                for (int i = 0; i < perRankDataSend[r]; ++i)
                    send_buf_data[r][i] = buf_x[send_col_idx[r][i]];
                MPI_Isend(send_buf_data[r], perRankDataSend[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD,
                          &send_reqs[r]);
            } else
                send_reqs[r] = MPI_REQUEST_NULL;
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
            recvColFromRanks[r][reqColFromRank[r][i] - procs_info[r].first_row] = recv_buf[r][i];
        }
    }

    for (int i = 0; i < procs_info[rank].M; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            /// Local elements multiplication
            if (in_diagonal(col_ptr[k], procs_info[rank].first_row, procs_info[rank].last_row)) {
                y[i] += val_ptr[k] * buf_x[col_ptr[k] - procs_info[rank].first_row];
            } else {
                /// Global elements multiplication
                int r = getRank(nRanks, procs_info, col_ptr[k]);
                y[i] += val_ptr[k] * recvColFromRanks[r][col_ptr[k] - procs_info[r].first_row];
            }
        }
    }
    /*
    if (reqRequired > 0) {
        /// Global elements multiplication
        for (int i = 0; i < procs_info[rank].M; ++i) {
            for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
                if (!in_diagonal(col_ptr[k], procs_info[rank].first_row, procs_info[rank].last_row)) {
                    int r = getRank(nRanks, procs_info, col_ptr[k]);
                    y[i] += val_ptr[k] * recvColFromRanks[r][col_ptr[k] - procs_info[r].first_row];
                }
            }
        }
    }
*/
    if(nRanksExpectCol > 0) {
        /// Wait until send request delivered to through network.
        MPI_Waitall(nRanks, send_reqs, MPI_STATUS_IGNORE);
    }
    if (reqRequired > 0 || nRanksExpectCol > 0) {
        for (int r = 0; r < nRanks; r++) {
            if (perRankDataRecv[r] > 0) {
                free(recv_buf[r]);
                free(recvColFromRanks[r]);
            }
            if (perRankDataSend[r] > 0)
                free(send_buf_data[r]);
        }
    }

    if (reqRequired > 0) {
        free(recv_buf);
        free(recvColFromRanks);
    }
    if (reqRequired > 0) {
        free(send_buf_data);
    }

    return y;
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

int findInterRanksComm(int rank, int nRanks, proc_info_t *procs_info, int *col_ptr, int *perRankDataRecv,
                       int **reqColFromRank,
                       int *count_communication, int *interProcessCall) {
    (*count_communication) = 0;
    (*interProcessCall) = 0;
    /* build sending blocks to processors */
    Map *map = (Map *) malloc_or_exit(procs_info[rank].NZ * sizeof(Map));
    int dest, col, reqRequired = 0;
    for (int i = 0; i < procs_info[rank].NZ; i++) {
        col = col_ptr[i];
        /// Check off-diagonal nonzero elements that belongs to other ranks
        if (in_diagonal(col, procs_info[rank].first_row, procs_info[rank].last_row) || !(getVal(map, col, i + 1) < 0))
            continue;
        ///search which rank has the element
        dest = -1;
        for (int r = 0; r < nRanks; r++) {
            if (in_diagonal(col, procs_info[r].first_row, procs_info[r].last_row)) {
                dest = r;
                break;
            }
        }
        assert(dest >= 0);
        reqRequired++;
        ///insert new request
        reqColFromRank[dest][perRankDataRecv[dest]++] = col;
        map[i].key.col = col;
        map[i].value.val = 1;
        if (perRankDataRecv[dest] == 1)
            (*interProcessCall)++;
        (*count_communication)++;
    }
    free(map);
    return reqRequired;
}

/**
 *
 * @param argc
 * @param argv
 * @return
 */
int shareReqColumnInfo(int rank, int nRanks, proc_info_t *procs_info, int *perRankDataRecv, int **reqColFromRank,
                       int *perRankDataSend, int **send_col_idx, int reqRequired) {
    /// Send Requests
    int *expect = (int *) calloc_or_exit(nRanks, sizeof(int));
    MPI_Request *send_reqs;
    if (reqRequired > 0) {
        send_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
        for (int r = 0; r < nRanks; r++) {
            if (r == rank || perRankDataRecv[r] == 0) {
                send_reqs[r] = MPI_REQUEST_NULL;
                continue;
            }
            expect[r] = 1;
            /// send the request
            MPI_Isend(reqColFromRank[r], perRankDataRecv[r], MPI_INT, r, REQUEST_TAG, MPI_COMM_WORLD, &send_reqs[r]);
        }
    }
    int *all_process_expect = (int *) calloc_or_exit(nRanks, sizeof(int));
    MPI_Allreduce(expect, all_process_expect, nRanks, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    /// Receive the requests
    int *reqs;
    MPI_Status status;
    int req_count;
    for (int p = 0; p < all_process_expect[rank]; p++) {
        /// Wait until a request comes
        MPI_Probe(MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &req_count);
        /// reply to this proccess
        int r = status.MPI_SOURCE;
        send_col_idx[r] = (int *) malloc_or_exit(req_count * sizeof(int));
        reqs = (int *) malloc_or_exit(req_count * sizeof(int));
        perRankDataSend[r] = req_count;
        MPI_Recv(reqs, req_count, MPI_INT, r, REQUEST_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < req_count; i++) {
            if (reqs[i] < procs_info[rank].first_row || reqs[i] > procs_info[rank].last_row) {
                printf("Wrong index %d looking at process %d\n", reqs[i], p);
                return 0;
            }
            send_col_idx[r][i] = reqs[i];
        }
    }
    if (reqRequired > 0)
        MPI_Waitall(nRanks, send_reqs, MPI_STATUS_IGNORE);
    if (all_process_expect[rank]>0 && reqs != NULL)
        free(reqs);
    if (expect != NULL)
        free(expect);
    return all_process_expect[rank];
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

    double comp_time = 0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0, sparsity = 0;;
    int nonZeroPerRow = 0, total_run = 100, mat_row = 0, mat_col = 0;
    int nRanks, rank;
    int *row_ptr, *col_ptr;
    double *val_ptr, *buf_x, *res;
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

    if (argc < 3) {
        printf("Usage: %s matrix_row matrix_col nonZeroPerRow [Total_Runs]\n", argv[0]);
        return 0;
    } else {
        mat_row = atoi(argv[1]);
        mat_col = atoi(argv[2]);
        nonZeroPerRow = atoi(argv[3]);
        if (argc > 4)
            total_run = atoi(argv[4]);
        if (argc > 5)
            sparsity = atoi(argv[5]);
    }
    ranks_info[rank].M = mat_row;
    ranks_info[rank].N = mat_col;
    ranks_info[rank].NZ = nonZeroPerRow * mat_row;;
    ranks_info[rank].first_row = rank * mat_row;
    ranks_info[rank].last_row = (rank + 1) * mat_row;
    if (nonZeroPerRow <= 0) {
        printf("[%d], There will must one non zero column in the matrix in every row\n", rank);
        return 0;
    }
    if (nonZeroPerRow > mat_row) {
        if (rank == MASTER) {
            printf("[%d] nonzero=%d, max nonzero=%d, number process=%d\n", rank, nonZeroPerRow, mat_row, nRanks);
        }
        nonZeroPerRow = mat_row;
        ranks_info[rank].NZ = nonZeroPerRow * mat_row;
    }
    if (ranks_info[rank].NZ <= 0) {
        printf("[%d] Matrix can not be sized zero=%d\n", rank, ranks_info[rank].NZ);
        return 0;
    }
    /// Initialize CSR row, col and value pointer.
    row_ptr = (int *) malloc((mat_row + 1) * sizeof(int));
    col_ptr = (int *) malloc(ranks_info[rank].NZ * sizeof(int));
    val_ptr = (double *) malloc(ranks_info[rank].NZ * sizeof(double));
    /// Create random CSR matrix with the given parameter
    if (csr_random_mat(rank, ranks_info, row_ptr, col_ptr, val_ptr, mat_row, mat_col, nonZeroPerRow, sparsity) != 1) {
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nRanks, (ranks_info[rank].M*nRanks),
               nonZeroPerRow);
    }

    /// Create CSR Diagonal matrix with the given parameter
    /*if (csr_diagonal_mat(rank, row_ptr, col_ptr, val_ptr, mat_row, nonZeroPerRow) != 1) {
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nRanks,
               (ranks_info[rank].M * nRanks),
               nonZeroPerRow);
    }*/

    /// Create vector x and fill with 1.0
    buf_x = (double *) malloc_or_exit(mat_row * sizeof(double));
    for (int i = 0; i < mat_row; i++) {
        buf_x[i] = 1.00;
    }
    /// Share process info among all the processes
    MPI_Allgather(&ranks_info[rank], 1, procs_info_type, procs_info, 1, procs_info_type, MPI_COMM_WORLD);
    int *perRankDataRecv, **reqColFromRank;
    perRankDataRecv = (int *) calloc_or_exit(nRanks, sizeof(int));
    /// Allocate buffers for requests sending
    reqColFromRank = (int **) malloc_or_exit(nRanks * sizeof(int *));
    for (int i = 0; i < nRanks; i++) {
        if (i != rank && procs_info[i].M > 0)
            reqColFromRank[i] = (int *) malloc_or_exit(procs_info[i].M * sizeof(int));
    }

    int count_communication = 0, interProcessCall = 0, totalInterProcessCall = 0, avg_communication = 0, per_rank_data_send = 0;
    /// Find the columns that belong to other ranks
    int reqRequired = findInterRanksComm(rank, nRanks, procs_info, col_ptr, perRankDataRecv, reqColFromRank,
                                         &count_communication, &interProcessCall);
    if (reqRequired > 0) {
        if (interProcessCall > 0)
            avg_communication = count_communication / interProcessCall;
        MPI_Reduce(&interProcessCall, &totalInterProcessCall, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
        MPI_Reduce(&avg_communication, &per_rank_data_send, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    }
    int *perRankDataSend = (int *) calloc_or_exit(nRanks, sizeof(int));
    int **send_col_idx = (int **) malloc_or_exit(nRanks * sizeof(int *));

    int nRanksExpectCol = shareReqColumnInfo(rank, nRanks, procs_info, perRankDataRecv, reqColFromRank, perRankDataSend,
                                             send_col_idx, reqRequired);

    /// Start sparse matrix vector multiplication for each rank
    MPI_Barrier(MPI_COMM_WORLD);
    double start_time = MPI_Wtime();
    for (int r = 0; r < total_run; ++r) {
        res = matMull(rank, procs_info, nRanks, row_ptr, col_ptr, val_ptr, buf_x, send_col_idx, perRankDataRecv,
                      reqColFromRank, perRankDataSend, reqRequired, nRanksExpectCol);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    comp_time = (MPI_Wtime() - start_time) * 1000.00;
    avg_time = comp_time / total_run;

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nRanks;

//    char *outputFIle = (char *) malloc_or_exit(100 * sizeof(char));
//    strcpy(outputFIle, "CSR_SpMV_Model.csv");
//    strcpy(outputFIle, "CSR_SpMV_Model_Diagonal.csv");
    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n", rank,
               min_time, max_time, mean, procs_info[rank].NZ);
        FILE *resultCSV;
        FILE *checkFile;
        if ((checkFile = fopen("CSR_SpMV_Model_Random.csv", "r")) != NULL) {
            // file exists
            fclose(checkFile);
            if (!(resultCSV = fopen("CSR_SpMV_Model_Random.csv", "a"))) {
                fprintf(stderr, "fopen: failed to open file CSR_SpMV_Model_Random.csv");
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen("CSR_SpMV_Model_Random.csv", "w"))) {
                fprintf(stderr, "fopen: failed to open file CSR_SpMV_Model_Random.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "MatrixSize,MinTime,MaxTime,AvgTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock,Sparsity,AvgCommunication,AvgInterProcessCall\n");
        }

        fprintf(resultCSV, "%d,%10.3lf,%10.3lf,%10.3lf,%d,%d,%d,%d,%d,%d,%d\n", procs_info[rank].N, min_time, max_time,
                mean,
                total_run, nRanks, nonZeroPerRow, procs_info[rank].NZ, sparsity, (per_rank_data_send / nRanks),
                (totalInterProcessCall / nRanks));
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file CSR_SpMV_Model_Random.csv");
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


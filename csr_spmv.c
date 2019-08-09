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

double *
matMull(int rank, proc_info_t *procs_info, int nRanks, int *row_ptr, int *col_ptr, double *val_ptr,
        int *off_row_ptr, int *off_col_ptr, double *off_val_ptr,
        double *buf_x, int **send_col_idx, int *perRankDataRecv, int **reqColFromRank,
        int *perRankDataSend, int reqRequired, int nRanksExpectCol) {

    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(procs_info[rank].M, sizeof(double));
/// receiving blocks storage
    double **recv_buf, *recvColFromRanks;
    MPI_Request *send_reqs, *recv_reqs;
    int reqMade = 0;
    if (reqRequired > 0) {
        recv_buf = (double **) malloc_or_exit(nRanks * sizeof(double *));
        recvColFromRanks = (double *) malloc_or_exit(procs_info[rank].N * sizeof(double));
        /// MPI request storage
        recv_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
        for (int r = 0; r < nRanks; ++r) {
            if (r == rank || perRankDataRecv[r] <= 0) {
                recv_reqs[r] = MPI_REQUEST_NULL;
                continue;
            }
            recv_buf[r] = (double *) malloc_or_exit(perRankDataRecv[r] * sizeof(double));
            /// Receive the block (when it comes)
            MPI_Irecv(recv_buf[r], perRankDataRecv[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &recv_reqs[r]);
            reqMade++;
        }
    }
    /// Local elements multiplication
    for (int i = 0; i < procs_info[rank].M; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k)
            y[i] += val_ptr[k] * buf_x[col_ptr[k] - procs_info[rank].first_row];
    }

    double **send_buf_data;
    if (nRanksExpectCol > 0) {
        send_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
        /// Reply to the requests.
        send_buf_data = (double **) malloc_or_exit(nRanks * sizeof(double *));
        for (int r = 0; r < nRanks; ++r) {
            if (r == rank || perRankDataSend[r] <= 0){
                send_reqs[r] = MPI_REQUEST_NULL;
                continue;
            }
            send_buf_data[r] = (double *) malloc_or_exit(perRankDataSend[r] * sizeof(double));
            if (send_col_idx[r] == NULL){
                printf("[%d] Sending column not found for=%d\n", rank, r);
            }
            for (int i = 0; i < perRankDataSend[r]; ++i) {
                if (send_col_idx[r][i] < procs_info[rank].first_row || send_col_idx[r][i] > procs_info[rank].last_row) {
                    printf("Wrong index %d looking at process %d\n", send_col_idx[r][i], r);
                    return 0;
                }
                send_buf_data[r][i] = buf_x[send_col_idx[r][i] - procs_info[r].first_row];
            }
            MPI_Isend(send_buf_data[r], perRankDataSend[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &send_reqs[r]);
        }
    }

    printf("[%d] Request send done\n", rank);
    return y;
    int r;
    for (int q = 0; q < reqMade; q++) {
        MPI_Waitany(nRanks, recv_reqs, &r, MPI_STATUS_IGNORE);
        assert(r != MPI_UNDEFINED);
        /// fill x array with new elements.
        for (int i = 0; i < perRankDataRecv[r]; i++) {
            if (reqColFromRank[r][i]<0 || reqColFromRank[r][i] >= procs_info[rank].N){
                printf("[%d] Column=%d out of range\n", rank, reqColFromRank[r][i]);
                return 0;
            }
            recvColFromRanks[reqColFromRank[r][i]] = recv_buf[r][i];
        }
    }
    printf("[%d] Global data received\n", rank);

    if (reqRequired > 0) {
        /// Global elements multiplication
        for (int i = 0; i < procs_info[rank].M; ++i) {
            for (int k = off_row_ptr[i]; k < off_row_ptr[i + 1]; ++k)
                y[i] += off_val_ptr[k] * recvColFromRanks[off_col_ptr[k]];
        }
    }
    printf("[%d] Global multiplication done\n", rank);
    if (nRanksExpectCol > 0) {
        /// Wait until send request delivered to through network.
        MPI_Waitall(nRanks, send_reqs, MPI_STATUS_IGNORE);
    }
    if (reqRequired > 0 || nRanksExpectCol > 0) {
        for (int r = 0; r < nRanks; r++) {
            if (perRankDataRecv[r] > 0) {
                free(recv_buf[r]);
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

int findInterRanksComm(int rank, int nRanks, proc_info_t *procs_info, int *col_ptr, int offDiagonalElements, int *perRankDataRecv,
                       int **reqColFromRank, int *count_communication, int *interProcessCall) {
    (*count_communication) = 0;
    (*interProcessCall) = 0;
    /* build sending blocks to processors */
    int dest, col, reqRequired = 0;
    int *map = (int *) calloc_or_exit(procs_info[rank].N, sizeof(int));
    for (int i = 0; i < offDiagonalElements; i++) {
        col = col_ptr[i];
        /// Check is already calculated
        if (map[col] > 0)
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
        map[col] = 1;
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
    if (all_process_expect[rank] > 0 && reqs != NULL)
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

    double comp_time = 0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    int nonZeroPerRow = 0, total_run = 100, mat_row = 0, mat_col = 0, sparsity = 0;
    int nRanks, rank;
    int *row_ptr, *on_diagonal_row, *off_diagonal_row, *col_ptr, *off_diagonal_col, *on_diagonal_col;
    double *val_ptr, *off_diagonal_val, *on_diagonal_val, *buf_x, *res;
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
    ranks_info[rank].last_row = (rank + 1) * mat_row-1;
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
    int offDiagonalElements = 0;
    /// Create random CSR matrix with the given parameter
    if (csr_random_mat(rank, row_ptr, col_ptr, val_ptr, mat_row, mat_col, nonZeroPerRow, sparsity, &offDiagonalElements) != 1) {
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nRanks, (ranks_info[rank].M*nRanks),
               nonZeroPerRow);
    }
    /*/// Create CSR Diagonal matrix with the given parameter
    if (csr_diagonal_mat(rank, row_ptr, col_ptr, val_ptr, mat_row, nonZeroPerRow) != 1) {
        printf("[%d] Matrix Creation Failed process=%d, matrix size=%d, nonzero=%d\n", rank, nRanks,
               (ranks_info[rank].M * nRanks),
               nonZeroPerRow);
    }*/

    int diagonal_elements = ranks_info[rank].NZ - offDiagonalElements;
    if(diagonal_elements>0) {
        on_diagonal_row = (int *) malloc_or_exit((mat_row + 1) * sizeof(int));
        on_diagonal_col = (int *) malloc_or_exit(diagonal_elements * sizeof(int));
        on_diagonal_val = (double *) malloc_or_exit(diagonal_elements * sizeof(double));
        on_diagonal_row[0] = 0;
    }
    if(offDiagonalElements>0) {
        off_diagonal_row = (int *) malloc_or_exit((mat_row + 1) * sizeof(int));
        off_diagonal_col = (int *) malloc_or_exit(offDiagonalElements * sizeof(int));
        off_diagonal_val = (double *) malloc_or_exit(offDiagonalElements * sizeof(double));
        off_diagonal_row[0] = 0;
    }
    int on_diag_idx = 0, off_diag_idx = 0;
    for (int k = 0; k < mat_row; ++k) {
        for (int l = row_ptr[k]; l < row_ptr[k + 1]; ++l) {
            if (in_diagonal(col_ptr[l], ranks_info[rank].first_row, ranks_info[rank].last_row)) {
                on_diagonal_col[on_diag_idx] = col_ptr[l];
                on_diagonal_val[on_diag_idx] = val_ptr[l];
                on_diag_idx++;
            } else {
                off_diagonal_col[off_diag_idx] = col_ptr[l];
                off_diagonal_val[off_diag_idx] = val_ptr[l];
                off_diag_idx++;
            }
        }
        if(diagonal_elements>0)
            on_diagonal_row[k + 1] = on_diag_idx;

        if(offDiagonalElements>0)
            off_diagonal_row[k + 1] = off_diag_idx;
    }

    free(row_ptr);
    free(col_ptr);
    free(val_ptr);
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

    int count_communication = 0, interProcessCall = 0, totalInterProcessCall = 0, avg_communication = 0, per_rank_data_send = 0, reqRequired = 0;
    /// Find the columns that belong to other ranks
    if (offDiagonalElements > 0)
        reqRequired = findInterRanksComm(rank, nRanks, procs_info, off_diagonal_col, offDiagonalElements, perRankDataRecv, reqColFromRank,
                                         &count_communication, &interProcessCall);
    for(int r=0; r<nRanks; ++r) {
        for (int i = 0; i < perRankDataRecv[r]; i++) {
            if (reqColFromRank[r][i] < 0 || reqColFromRank[r][i] >= procs_info[rank].N) {
                printf("[%d] Column=%d out of range\n", rank, reqColFromRank[r][i]);
                return 0;
            }
        }
    }
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


    /**
     * start
     */
//    MPI_Request *send_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
    /// Reply to the requests.
    double **send_buf_data = (double **) malloc_or_exit(nRanks * sizeof(double *));
    for (int r = 0; r < nRanks; ++r) {
        if (r == rank || perRankDataSend[r] <= 0){
//            send_reqs[r] = MPI_REQUEST_NULL;
            continue;
        }
        send_buf_data[r] = (double *) malloc_or_exit(perRankDataSend[r] * sizeof(double));
        if (send_col_idx[r] == NULL){
            printf("[%d] Sending column not found for=%d\n", rank, r);
        }
        for (int i = 0; i < perRankDataSend[r]; ++i) {
            if (send_col_idx[r][i] < procs_info[rank].first_row || send_col_idx[r][i] > procs_info[rank].last_row) {
                printf("Wrong index %d looking at process %d\n", send_col_idx[r][i], r);
                return 0;
            }
            double col_val = buf_x[send_col_idx[r][i] - procs_info[r].first_row];
            send_buf_data[r][i] = (double) 120.0;
            printf("[%d] col value=%lf", col_val);
//            send_buf_data[r][i] = buf_x[send_col_idx[r][i] - procs_info[r].first_row];
        }
//        MPI_Isend(send_buf_data[r], perRankDataSend[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &send_reqs[r]);
    }
    printf("[%d] Done!! no problem...", rank);
    /**
     * end
     */
    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
    /// Start sparse matrix vector multiplication for each rank
    double start_time = MPI_Wtime();
    for (int r = 0; r < total_run; ++r) {
        res = matMull(rank, procs_info, nRanks, on_diagonal_row, on_diagonal_col, on_diagonal_val, off_diagonal_row,
                      off_diagonal_col, off_diagonal_val, buf_x, send_col_idx, perRankDataRecv, reqColFromRank,
                      perRankDataSend, reqRequired, nRanksExpectCol);
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
//    strcpy(outputFIle, "CSR_SpMV_Model_Random.csv");
    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d, Sparsity: %d\n",
               rank,
               min_time, max_time, mean, procs_info[rank].NZ, sparsity);
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
    free(on_diagonal_row);
    free(on_diagonal_col);
    free(on_diagonal_val);
    if(offDiagonalElements>0) {
        free(off_diagonal_row);
        free(off_diagonal_col);
        free(off_diagonal_val);
    }
    free(buf_x);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


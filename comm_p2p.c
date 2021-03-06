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

double *
communication(int rank, proc_info_t *procs_info, int nRanks, double *off_val_ptr, double *buf_x, int **send_col_idx,
              int *perRankDataRecv, int *colCount, int **reqColFromRank, int ***reqRowCol, int *perRankDataSend,
              int reqRequired, int nRanksExpectCol, double **recv_buf, double **send_buf_data, MPI_Request *recv_reqs,
              MPI_Request *send_reqs) {

    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(procs_info[rank].M, sizeof(double));
/// receiving blocks storage
    int reqMade = 0;
    if (reqRequired > 0) {
        for (int r = 0; r < nRanks; ++r) {
            if (r == rank || perRankDataRecv[r] <= 0) {
                continue;
            }
            /// Receive the block (when it comes)
            MPI_Irecv(recv_buf[r], perRankDataRecv[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &recv_reqs[r]);
            reqMade++;
        }
    }

    if (nRanksExpectCol > 0) {
        for (int r = 0; r < nRanks; ++r) {
            if (r == rank || perRankDataSend[r] <= 0) {
                continue;
            }
            for (int i = 0; i < perRankDataSend[r]; ++i) {
                send_buf_data[r][i] = buf_x[send_col_idx[r][i] - procs_info[rank].first_row];
            }
            MPI_Isend(send_buf_data[r], perRankDataSend[r], MPI_DOUBLE, r, RECEIVE_TAG, MPI_COMM_WORLD, &send_reqs[r]);
        }
    }

    int r;
    for (int q = 0; q < reqMade; q++) {
        MPI_Waitany(nRanks, recv_reqs, &r, MPI_STATUS_IGNORE);
        assert(r != MPI_UNDEFINED);
        /// fill x array with new elements.
        for (int i = 0; i < perRankDataRecv[r]; i++) {
            for (int j = 1; j < (1 + (2 * colCount[reqColFromRank[r][i]])); j += 2) {
                y[reqRowCol[r][i][j]] += off_val_ptr[reqRowCol[r][i][j + 1]] * recv_buf[r][i];
            }
        }
    }

    if (nRanksExpectCol > 0) {
        /// Wait until send request delivered to through network.
        MPI_Waitall(nRanks, send_reqs, MPI_STATUS_IGNORE);
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

int
getRemoteColumnInfo(int rank, int nRanks, proc_info_t *procs_info, int *row_ptr, int *col_ptr, int offDiagonalElements,
                    int *perRankDataRecv, int *colCount,
                    int **reqColFromRank, int ***reqRowCol, int *count_communication, int *interProcessCall) {
    (*count_communication) = 0;
    (*interProcessCall) = 0;
    /// build sending blocks to processors
    int dest, col, reqRequired = 0;
    int *map = (int *) calloc_or_exit(procs_info[rank].N, sizeof(int));
    int *trackColIdx = (int *) malloc_or_exit(procs_info[rank].N * sizeof(int));
    int *colWiseRank = (int *) calloc_or_exit(procs_info[rank].N, sizeof(int));
    int *perRankColCount = (int *) calloc_or_exit(procs_info[rank].N, sizeof(int));
    int *perColCount = (int *) calloc_or_exit(procs_info[rank].N, sizeof(int));
    for (int i = 0; i < offDiagonalElements; i++) {
        col = col_ptr[i];
        colCount[col] += 1;
        trackColIdx[col] = -1;
        if (map[col] > 0) {
            continue;
        }
        ///search which rank has the element
        dest = -1;
        for (int r = 0; r < nRanks; r++) {
            if (in_diagonal(col, procs_info[r].first_row, procs_info[r].last_row)) {
                dest = r;
                break;
            }
        }
        assert(dest >= 0);
        colWiseRank[col] = dest;
        reqRequired++;
        ///insert new request
        reqColFromRank[dest][perRankDataRecv[dest]++] = col;
        map[col] = 1;
        if (perRankDataRecv[dest] == 1)
            (*interProcessCall)++;
        (*count_communication)++;
    }
    for (int r = 0; r < nRanks; ++r) {
        reqRowCol[r] = (int **) malloc_or_exit(perRankDataRecv[r] * sizeof(int *));
        for (int i = 0; i < perRankDataRecv[r]; ++i) {
            reqRowCol[r][i] = (int *) malloc_or_exit((1 + (2 * colCount[reqColFromRank[r][i]])) * sizeof(int));
        }
    }

    for (int i = 0; i < procs_info[rank].M; ++i) {
        for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
            col = col_ptr[k];
            dest = colWiseRank[col];
            if (trackColIdx[col] < 0) {
                trackColIdx[col] = perRankColCount[dest]++;
            }
            if (perColCount[col] == 0) {
                reqRowCol[dest][trackColIdx[col]][perColCount[col]++] = col;
            }
            reqRowCol[dest][trackColIdx[col]][perColCount[col]++] = i;
            reqRowCol[dest][trackColIdx[col]][perColCount[col]++] = k;
        }
    }
    for (int r = 0; r < nRanks; ++r) {
        if (r == rank)
            continue;
        for (int i = 0; i < perRankDataRecv[r]; ++i) {
            if (reqColFromRank[r][i] != reqRowCol[r][i][0]) {
                printf("[%d] Data %d is not same to %d\n", r, reqColFromRank[r][i], reqRowCol[r][i][0]);
            }
        }
    }

    free(map);
    free(colWiseRank);
    free(perRankColCount);
    free(perColCount);
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

    char *in_file, *out_file = NULL;
    double comp_time = 0, min_time = 0.0, max_time = 0.0, avg_time = 0.0, mean = 0.0;
    int total_run = 100, skip=100, nRanks, rank, knl=0;
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
    if (diagonal_elements > 0) {
        on_diagonal_row = (int *) malloc_or_exit((ranks_info[rank].M + 1) * sizeof(int));
        on_diagonal_col = (int *) malloc_or_exit(diagonal_elements * sizeof(int));
        on_diagonal_val = (double *) malloc_or_exit(diagonal_elements * sizeof(double));
        on_diagonal_row[0] = 0;
    }
    if (offDiagonalElements > 0) {
        off_diagonal_row = (int *) malloc_or_exit((ranks_info[rank].M + 1) * sizeof(int));
        off_diagonal_col = (int *) malloc_or_exit(offDiagonalElements * sizeof(int));
        off_diagonal_val = (double *) malloc_or_exit(offDiagonalElements * sizeof(double));
        off_diagonal_row[0] = 0;
    }
    int on_diag_idx = 0, off_diag_idx = 0;
    int in_diagonal_bandwidth = 0, bandwidth = 0;
    for (int k = 0; k < ranks_info[rank].M; ++k) {
        int l_col = ranks_info[rank].last_row, h_col = ranks_info[rank].first_row, max_col =ranks_info[rank].last_row, min_col=ranks_info[rank].first_row;
        for (int l = row_ptr[k]; l < row_ptr[k + 1]; ++l) {
            if (in_diagonal(col_ptr[l], ranks_info[rank].first_row, ranks_info[rank].last_row)) {
                on_diagonal_col[on_diag_idx] = col_ptr[l];
                on_diagonal_val[on_diag_idx] = val_ptr[l];
                on_diag_idx++;
                if (col_ptr[l] < l_col) {
                    l_col = col_ptr[l];
                }
                if (col_ptr[l] > h_col) {
                    h_col = col_ptr[l];
                }
            } else {
                off_diagonal_col[off_diag_idx] = col_ptr[l];
                off_diagonal_val[off_diag_idx] = val_ptr[l];
                off_diag_idx++;
            }
            if(col_ptr[l]<min_col){
                min_col = col_ptr[l];
            }
            if(col_ptr[l]>max_col){
                max_col = col_ptr[l];
            }
        }
        if (in_diagonal_bandwidth < (h_col - l_col + 1)) {
            in_diagonal_bandwidth = (h_col - l_col + 1);
        }
        if(bandwidth < (max_col - min_col + 1)){
            bandwidth = (max_col - min_col + 1);
        }
        if (diagonal_elements > 0)
            on_diagonal_row[k + 1] = on_diag_idx;

        if (offDiagonalElements > 0)
            off_diagonal_row[k + 1] = off_diag_idx;
    }
    free(row_ptr);
    free(col_ptr);
    free(val_ptr);
    /// Create vector x and fill with 1.0
    buf_x = (double *) malloc_or_exit(ranks_info[rank].M * sizeof(double));
    for (int i = 0; i < ranks_info[rank].M; i++) {
        buf_x[i] = 1.00;
    }
    /// Share process info among all the processes
    MPI_Allgather(&ranks_info[rank], 1, procs_info_type, procs_info, 1, procs_info_type, MPI_COMM_WORLD);
    int *perRankDataRecv, **reqColFromRank, ***reqRowCol;
    int *colCount = (int *) calloc_or_exit(procs_info[rank].N, sizeof(int));
    perRankDataRecv = (int *) calloc_or_exit(nRanks, sizeof(int));
    /// Allocate buffers for requests sending
    reqColFromRank = (int **) malloc_or_exit(nRanks * sizeof(int *));
    reqRowCol = (int ***) malloc_or_exit(nRanks * sizeof(int **));
    for (int i = 0; i < nRanks; i++) {
        if (i != rank && procs_info[i].M > 0) {
            reqColFromRank[i] = (int *) malloc_or_exit(procs_info[i].M * sizeof(int));
        }
    }

    int count_communication = 0, interProcessCall = 0, totalInterProcessCall = 0, avg_communication = 0, per_rank_data_send = 0, reqRequired = 0;
    /// Find the columns that belong to other ranks
    if (offDiagonalElements > 0) {
        /*reqRequired = findInterRanksComm(rank, nRanks, procs_info, off_diagonal_col, offDiagonalElements,
                                         perRankDataRecv, reqColFromRank,
                                         &count_communication, &interProcessCall);*/

        reqRequired = getRemoteColumnInfo(rank, nRanks, procs_info, off_diagonal_row, off_diagonal_col,
                                          offDiagonalElements,
                                          perRankDataRecv, colCount, reqColFromRank, reqRowCol, &count_communication,
                                          &interProcessCall);
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

    double **recv_buf, **send_buf_data;
    MPI_Request *send_reqs, *recv_reqs;
    if (reqRequired > 0) {
        recv_buf = (double **) malloc_or_exit(nRanks * sizeof(double *));
        /// MPI request storage
        recv_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
        for (int r = 0; r < nRanks; ++r) {
            if (r == rank || perRankDataRecv[r] <= 0) {
                recv_reqs[r] = MPI_REQUEST_NULL;
                continue;
            }
            recv_buf[r] = (double *) malloc_or_exit(perRankDataRecv[r] * sizeof(double));
        }
    }
    if (nRanksExpectCol > 0) {
        send_reqs = (MPI_Request *) malloc_or_exit(nRanks * sizeof(MPI_Request));
        /// Reply to the requests.
        send_buf_data = (double **) malloc_or_exit(nRanks * sizeof(double *));
        for (int r = 0; r < nRanks; ++r) {
            if (r == rank || perRankDataSend[r] <= 0) {
                send_reqs[r] = MPI_REQUEST_NULL;
                continue;
            }
            send_buf_data[r] = (double *) malloc_or_exit(perRankDataSend[r] * sizeof(double));
        }
    }

    /// Start sparse matrix vector multiplication for each rank
    double start_time = 0.0;
    struct timespec start, end;
    MPI_Barrier(MPI_COMM_WORLD);
    for (int r = 0; r < total_run; ++r) {
        clock_gettime(CLOCK_MONOTONIC, &start);
        res = communication(rank, procs_info, nRanks, off_diagonal_val, buf_x, send_col_idx, perRankDataRecv, colCount,
                            reqColFromRank, reqRowCol, perRankDataSend, reqRequired, nRanksExpectCol, recv_buf,
                            send_buf_data, recv_reqs, send_reqs);
        if(r>=skip) {
            clock_gettime(CLOCK_MONOTONIC, &end);
            comp_time += ((end.tv_sec * 1000 + (end.tv_nsec / 1.0e6)) - (start.tv_sec * 1000 + (start.tv_nsec / 1.0e6)));
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
//    comp_time = (MPI_Wtime() - start_time) * 1000.00;
    avg_time = comp_time / total_run;

    MPI_Reduce(&avg_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&avg_time, &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = mean / nRanks;

    int max_nnz = 0, max_row = 0, max_band_width = 0, max_in_diag_band_width = 0;
    double max_nnz_per_row = 0.0, nnz_per_row = procs_info[rank].NZ / procs_info[rank].M;
    MPI_Reduce(&procs_info[rank].NZ, &max_nnz, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&procs_info[rank].M, &max_row, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&in_diagonal_bandwidth, &max_in_diag_band_width, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&bandwidth, &max_band_width, 1, MPI_INT, MPI_MAX, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&nnz_per_row, &max_nnz_per_row, 1, MPI_DOUBLE, MPI_MAX, MASTER, MPI_COMM_WORLD);

    /// print execution stats
    if (rank == MASTER) {
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms, NonZero: %d\n",
               rank, min_time, max_time, mean, procs_info[rank].NZ);
        char *_ptr = strtok(in_file, "/");
        char *matrixName = strtok(strtok(NULL, "-"), ".");
        char outputFile[100] = "Skylake_CSR_Comm_on_MPI.csv";
        if(knl > 0)
            strcpy(outputFile, "KNL_CSR_Comm_on_MPI.csv");

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
                    "Name,MatrixSize,MaxRow,MinTime,MaxTime,AvgTime,TotalRun,nProcess,NonZeroPerRow,NonZeroPerBlock,AvgCommunication,AvgInterProcessCall,SizeOfData,DiagonalElem,OffDiagonalElem,Bandwidth,InDiagonalBandwidth\n");
        }

        fprintf(resultCSV, "%s,%d,%d,%lf,%lf,%lf,%d,%d,%lf,%d,%lf,%lf,%d,%d,%d,%d,%d\n", matrixName,
                procs_info[rank].N, max_row, min_time, max_time,
                mean, total_run, nRanks, max_nnz_per_row, max_nnz, ((double) per_rank_data_send / nRanks),
                ((double) totalInterProcessCall / nRanks), sizeof(double), diagonal_elements, offDiagonalElements, max_band_width, max_in_diag_band_width);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file %s", outputFile);
            exit(EXIT_FAILURE);
        }
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
    }
    if (reqRequired > 0) {
        free(send_buf_data);
    }

    if (diagonal_elements > 0) {
        free(on_diagonal_row);
        free(on_diagonal_col);
        free(on_diagonal_val);
    }
    if (offDiagonalElements > 0) {
        free(off_diagonal_row);
        free(off_diagonal_col);
        free(off_diagonal_val);
    }
    free(buf_x);
    free(colCount);
    /* MPI: end */
    MPI_Finalize();

    return 0;
}


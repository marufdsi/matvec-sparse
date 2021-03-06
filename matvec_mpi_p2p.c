#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

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
int *to_send, **send_buf;

enum tag {
    REQUEST_TAG, REPLY_TAG
};

double *mat_vec_mult_parallel(int rank, int nprocs, int *buf_i_idx, int *buf_j_idx, double *buf_values,
                              double *buf_x, int **rep_col_idx, int *expected_col, int req_made) {

    /// y vector for y = M*x/
    double *y = (double *) calloc_or_exit(proc_info[rank].M, sizeof(double));
    /// receiving blocks storage
    double **recv_buf = (double **) malloc_or_exit(nprocs * sizeof(double));
    for (int p = 0; p < nprocs; p++)
        if (to_send[p] > 0)
            recv_buf[p] = (double *) malloc_or_exit(to_send[p] * sizeof(double));

    /// MPI request storage
    MPI_Request *send_reqs = (MPI_Request *) malloc_or_exit(nprocs * sizeof(MPI_Request));
    MPI_Request *recv_reqs = (MPI_Request *) malloc_or_exit(nprocs * sizeof(MPI_Request));
    /// sending receive requests to processes in blocks
    for (int p = 0; p < nprocs; p++) {
        if (p == rank || to_send[p] == 0) {
            recv_reqs[p] = MPI_REQUEST_NULL;
            continue;
        }
        /// Receive the block (when it comes)
        MPI_Irecv(recv_buf[p], to_send[p], MPI_DOUBLE, p, REPLY_TAG, MPI_COMM_WORLD, &recv_reqs[p]);
    }

    /// Reply to the requests.
    double **rep_buf_data = (double **) malloc_or_exit(nprocs * sizeof(double));
    for (int p = 0; p < nprocs; ++p) {
        if (expected_col[p] > 0) {
            rep_buf_data[p] = (double *) malloc_or_exit(expected_col[p] * sizeof(double));
            for (int i = 0; i < expected_col[p]; ++i)
                rep_buf_data[p][i] = buf_x[rep_col_idx[p][i]];
            MPI_Isend(rep_buf_data[p], expected_col[p], MPI_DOUBLE, p, REPLY_TAG, MPI_COMM_WORLD, &send_reqs[p]);
        } else
            send_reqs[p] = MPI_REQUEST_NULL;
    }

    /// Local elements multiplication
    for (int k = 0; k < proc_info[rank].NZ; k++)
        if (in_diagonal(buf_j_idx[k], proc_info[rank].first_row, proc_info[rank].last_row))
            y[buf_i_idx[k] - proc_info[rank].first_row] +=
                    buf_values[k] * buf_x[buf_j_idx[k] - proc_info[rank].first_row];

    int p;
    /// need to update the initialization
    double *vecFromRemotePros = (double *) calloc_or_exit(proc_info[rank].N, sizeof(double));
    /// Wait to receive the request of the global column.
    for (int q = 0; q < req_made; q++) {
        MPI_Waitany(nprocs, recv_reqs, &p, MPI_STATUS_IGNORE);
        assert(p != MPI_UNDEFINED);

        /// fill x array with new elements.
        for (int i = 0; i < to_send[p]; i++)
            vecFromRemotePros[send_buf[p][i]] = recv_buf[p][i];
    }

    /// Global elements multiplication
    for (int k = 0; k < proc_info[rank].NZ; k++)
        if (!in_diagonal(buf_j_idx[k], proc_info[rank].first_row, proc_info[rank].last_row))
            y[buf_i_idx[k] - proc_info[rank].first_row] += buf_values[k] * vecFromRemotePros[buf_j_idx[k]];

    /// Wait until send request delivered to through network.
    MPI_Waitall(nprocs, send_reqs, MPI_STATUS_IGNORE);
    for (int p = 0; p < nprocs; p++) {
        if (to_send[p] > 0)
            free(recv_buf[p]);
        if (expected_col[p] > 0)
            free(rep_buf_data[p]);
    }

    free(recv_buf);
    free(rep_buf_data);
    free(vecFromRemotePros);
    return y;
}

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

int *CalculateInterProcessComm(int rank, int nprocs, int *buf_j_idx) {
    int count_communication = 0;
    int interProcessCall = 0;
    /* build sending blocks to processors */
    int *map = (int *) calloc_or_exit(proc_info[rank].N, sizeof(int));
    int dest, col;
    for (int i = 0; i < proc_info[rank].NZ; i++) {
        col = buf_j_idx[i];
        /// check whether I need to send a request
        if (in_diagonal(col, proc_info[rank].first_row, proc_info[rank].last_row) || map[col] > 0) {
            continue;
        }

        ///search which process has the element
        ////* NOTE: Due to small number or processes, serial search is faster

        dest = -1;
        for (int p = 0; p < nprocs; p++) {
            if (in_diagonal(col, proc_info[p].first_row, proc_info[p].last_row)) {
                dest = p;
                break;
            }
        }
        assert(dest >= 0);
        ///insert new request
        send_buf[dest][to_send[dest]++] = col;
        map[col] = 1;
        if (to_send[dest] == 1)
            interProcessCall++;
        count_communication++;
    }
    int total_communication, totalInterProcessCall, avg_communication = 0, per_rank_data_send;
    if (interProcessCall > 0)
        avg_communication = count_communication / interProcessCall;
    MPI_Reduce(&count_communication, &total_communication, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&interProcessCall, &totalInterProcessCall, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);
    MPI_Reduce(&avg_communication, &per_rank_data_send, 1, MPI_INT, MPI_SUM, MASTER, MPI_COMM_WORLD);

    free(map);
    int *returnPtr;
    if (rank == MASTER) {
        returnPtr = (int *) malloc_or_exit(3 * sizeof(int));
        returnPtr[0] = totalInterProcessCall / nprocs;
        returnPtr[1] = total_communication / nprocs;
        returnPtr[2] = per_rank_data_send / nprocs;
    }
    return returnPtr;
}

int main(int argc, char *argv[]) {
    char *in_file, *out_file = NULL;

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
    to_send = (int *) calloc_or_exit(nprocs, sizeof(int));
    /// Allocate buffers for requests sending
    send_buf = (int **) malloc_or_exit(nprocs * sizeof(int *));
    /// Read input matrix
    if (rank_wise_read_matrix(in_file, &buf_i_idx, &buf_j_idx, &buf_values,
                              &proc_info[rank].M, &proc_info[rank].N, &proc_info[rank].NZ,
                              &proc_info[rank].first_row, &proc_info[rank].last_row, rank) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }

    buf_x = (double *) malloc_or_exit(proc_info[rank].M * sizeof(double));
    for (int i = 0; i < proc_info[rank].M; i++) {
        buf_x[i] = 1;
    }

    /// Share process info among all the processes
    MPI_Allgather(&proc_info[rank], 1, proc_info_type, proc_info, 1, proc_info_type, MPI_COMM_WORLD);

    int *row_count, *row_offset;
    row_count = (int *) malloc_or_exit(nprocs * sizeof(int));
    row_offset = (int *) malloc_or_exit(nprocs * sizeof(int));
    for (int p = 0; p < nprocs; p++) {
        row_count[p] = proc_info[p].M;
        row_offset[p] = proc_info[p].first_row;
    }

    for (int i = 0; i < nprocs; i++) {
        if (i != rank && proc_info[i].M > 0)
            send_buf[i] = (int *) malloc_or_exit(proc_info[i].M * sizeof(int));
    }

    int *total_comm = CalculateInterProcessComm(rank, nprocs, buf_j_idx);

    /// Find the process that need data from the current rank.
    int *expect = (int *) calloc_or_exit(nprocs, sizeof(int));
    MPI_Request *send_reqs = (MPI_Request *) malloc_or_exit(nprocs * sizeof(MPI_Request));
    int req_made = 0;
    for (int p = 0; p < nprocs; p++) {
        /* need to send to this proc? */
        if (p == rank || to_send[p] == 0) {
            send_reqs[p] = MPI_REQUEST_NULL;
            continue;
        }
        /* logistics */
        req_made++;
        expect[p] = 1;
        /// send the request
        MPI_Isend(send_buf[p], to_send[p], MPI_INT, p, REQUEST_TAG, MPI_COMM_WORLD, &send_reqs[p]);
    }
    int *all_process_expect = (int *) calloc_or_exit(nprocs, sizeof(int));
    MPI_Allreduce(expect, all_process_expect, nprocs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if (rank == MASTER) {
        printf("[%d] Average inter rank MPI send requests=%d\n", rank, total_comm[0]);
        printf("[%d] Average Inter Rank Communication Required: %d\n", rank, total_comm[1]);
        printf("[%d] Average data send to a rank=%d\n", rank, total_comm[2]);
    }

    /// Receive the requests
    int *reqs;
    int *expected_col = (int *) calloc_or_exit(nprocs, sizeof(int));
    int **rep_col_idx = (int **) malloc_or_exit(nprocs * sizeof(int *)); /* reply blocks storage */
    MPI_Status status;
    int req_count;
    for (int p = 0; p < all_process_expect[rank]; p++) {
        /// Wait until a request comes
        MPI_Probe(MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &req_count);
        /// reply to this proccess
        int r_p = status.MPI_SOURCE;
        rep_col_idx[r_p] = (int *) malloc_or_exit(req_count * sizeof(int));
        reqs = (int *) malloc_or_exit(req_count * sizeof(int));
        expected_col[r_p] = req_count;
        MPI_Recv(reqs, req_count, MPI_INT, r_p, REQUEST_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < req_count; i++) {
            if (reqs[i] < proc_info[rank].first_row || reqs[i] > proc_info[rank].last_row) {
                printf("Wrong index %d looking at process %d\n", reqs[i], p);
                return 0;
            }
            rep_col_idx[r_p][i] = reqs[i] - proc_info[rank].first_row;
        }
    }
    MPI_Waitall(nprocs, send_reqs, MPI_STATUS_IGNORE);

    /* Matrix-vector multiplication for each processes */
    double min_time = 0, max_time, avg_time;
    double *y;
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    y = mat_vec_mult_parallel(rank, nprocs, buf_i_idx, buf_j_idx, buf_values, buf_x, rep_col_idx, expected_col,
                              req_made);
    double runTime = (MPI_Wtime() - t) * 1000.00;
    MPI_Barrier(MPI_COMM_WORLD);
    /// Calculate Execution Time
    MPI_Reduce(&runTime, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&runTime, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&runTime, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    avg_time = avg_time / nprocs;
    double compareTime = avg_time;
    if (rank == MASTER) {
        printf("[%d] First run MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms\n", rank, min_time, max_time,
               avg_time);
    }

    double totalTime = 0, mean = 0, latency = 0, timeRequired = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    t = MPI_Wtime();
    for (int r = 0; r < TOTAL_RUNS; r++) {
        y = mat_vec_mult_parallel(rank, nprocs, buf_i_idx, buf_j_idx, buf_values, buf_x, rep_col_idx, expected_col,
                                  req_made);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    totalTime = (MPI_Wtime() - t) * 1000.00;
    latency = totalTime / TOTAL_RUNS;
    MPI_Reduce(&latency, &min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&latency, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&latency, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    mean = avg_time / nprocs;
    if (rank == MASTER) {
        printf("[%d] Total Time: %10.3lf, Total Runs: %d\n", rank, totalTime, TOTAL_RUNS);
        printf("[%d] Min Latency: %10.3lf, Max Latency: %10.3lf, Avg Latency: %10.3lf\n", rank, min_time, max_time,
               mean);
    }

    int minNonZero = 0, maxNonZero = 0, avgNonZero;
    MPI_Reduce(&proc_info[rank].NZ, &minNonZero, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&proc_info[rank].NZ, &maxNonZero, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&proc_info[rank].NZ, &avgNonZero, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    avgNonZero = avgNonZero/nprocs;

    int minRequests = 0, maxRequests = 0, avgRequests = 0;
    MPI_Reduce(&all_process_expect[rank], &minRequests, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&all_process_expect[rank], &maxRequests, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&all_process_expect[rank], &avgRequests, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    avgRequests = avgRequests / nprocs;
    /// Print execution stats
    if (rank == MASTER) {
        char *_ptr = strtok(in_file, "_");
        printf("[%d] Computation MinTime: %10.3lf, MaxTime: %10.3lf, AvgTime: %10.3lf ms\n", rank, min_time, max_time,
               mean);
        FILE *resultCSV;
        FILE *checkFile;
        if ((checkFile = fopen("MPISpMVResult.csv", "r")) != NULL) {
            /// file exists
            fclose(checkFile);
            if (!(resultCSV = fopen("MPISpMVResult.csv", "a"))) {
                fprintf(stderr, "fopen: failed to open file MPISpMVResult.csv");
                exit(EXIT_FAILURE);
            }
        } else {
            if (!(resultCSV = fopen("MPISpMVResult.csv", "w"))) {
                fprintf(stderr, "fopen: failed to open file MPISpMVResult.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV,
                    "MatrixName,MinTime,MaxTime,AvgTime,TotalRun,nProcess,AvgSendRequests,AvgDataRequests,PerRankDataRequests,SizeOfMsg,MaxNonZero,MinNonZero,AvgNonZero,MinRequests,MaxRequests,AvgRequests\n");
        }
        fprintf(resultCSV, "%s,%10.3lf,%10.3lf,%10.3lf,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n", _ptr, min_time, max_time, mean,
                TOTAL_RUNS, nprocs, total_comm[0], total_comm[1], total_comm[2], sizeof(double), maxNonZero, minNonZero,avgNonZero,
                minRequests, maxRequests, avgRequests);
        if (fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file MPISpMVResult");
            exit(EXIT_FAILURE);
        }
    }

    free(buf_values);
    free(buf_i_idx);
    free(buf_j_idx);
    free(buf_x);
    free(rep_col_idx);
    free(expected_col);
    free(expect);
    free(send_buf);
    free(to_send);
    free(all_process_expect);

    /* MPI: end */
    MPI_Finalize();

    return 0;
}


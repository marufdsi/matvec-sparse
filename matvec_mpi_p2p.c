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
//enum policies policy = EQUAL_NZ;
enum policies policy = EQUAL_ROWS;
MPI_Datatype proc_info_type;
proc_info_t *proc_info;

enum tag {
    REQUEST_TAG, REPLY_TAG
};

double *mat_vec_mult_parallel(int rank, int nprocs, int *buf_i_idx, int *buf_j_idx, double *buf_values, double *buf_x,
                              int *row_count, int *row_offset/*, int **send_buf, int *to_send*/) {
    double *res;            /* result of multiplication res = A*x */

    /* allocate memory for vectors and submatrixes */
    double *y = (double *) calloc_or_exit(proc_info[rank].M, sizeof(double));
    if (rank == MASTER) {
        res = (double *) malloc_or_exit(proc_info[rank].N * sizeof(double));
    }

    /* allocate buffers for requests sending */
     int **send_buf = (int **) malloc_or_exit(nprocs * sizeof(int *));
     for (int i = 0; i < nprocs; i++) {
         if (i != rank && proc_info[i].M > 0)
             send_buf[i] = (int *) malloc_or_exit(proc_info[i].M * sizeof(int));
     }

     int *to_send = (int *) calloc_or_exit(nprocs, sizeof(int));    // # of req to each proc
    int *map = (int *) calloc_or_exit(proc_info[rank].N, sizeof(int));

    /// build sending blocks to processors
    int dest, col;
    for (int i = 0; i < proc_info[rank].NZ; i++) {
        col = buf_j_idx[i];
       ///  check whether I need to send a request
        if (in_diagonal(col, proc_info[rank].first_row, proc_info[rank].last_row) || map[col] > 0) {
            continue;
        }

         ///search which process has the element
         ///* NOTE: Due to small number or processes, serial search is faster

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
    }

    /* MPI request storage */
    MPI_Request *send_reqs = (MPI_Request *) malloc_or_exit(nprocs * sizeof(MPI_Request));
    MPI_Request *recv_reqs = (MPI_Request *) malloc_or_exit(nprocs * sizeof(MPI_Request));

    /* receiving blocks storage */
    double **recv_buf = (double **) malloc_or_exit(nprocs * sizeof(double *));
    for (int p = 0; p < nprocs; p++) {
        if (to_send[p] > 0)
            recv_buf[p] = (double *) malloc_or_exit(to_send[p] * sizeof(double));
    }

    /* sending requests to processes in blocks */
    int req_made = 0;
    int *expect = (int *) calloc_or_exit(nprocs, sizeof(int));
    for (int p = 0; p < nprocs; p++) {
        /* need to send to this proc? */
        if (p == rank || to_send[p] == 0) {
            send_reqs[p] = recv_reqs[p] = MPI_REQUEST_NULL;
            continue;
        }
        debug("[%d] Sending requests to process %2d \t[%5d]\n", rank, p, to_send[p]);

        /* logistics */
        expect[p] = 1;
        req_made++;

        /* send the request */
        MPI_Isend(send_buf[p], to_send[p], MPI_INT, p, REQUEST_TAG,
                  MPI_COMM_WORLD, &send_reqs[p]);
        /* recv the block (when it comes) */
        MPI_Irecv(recv_buf[p], to_send[p], MPI_DOUBLE, p, REPLY_TAG,
                  MPI_COMM_WORLD, &recv_reqs[p]);
    }

    int *all_process_expect = (int *) calloc_or_exit(nprocs, sizeof(int));
    MPI_Allreduce(expect, all_process_expect, nprocs, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    /**** reply to requests ****/
    int *reqs = (int *) malloc_or_exit(proc_info[rank].M * sizeof(int));
    double **rep_buf = (double **) malloc_or_exit(nprocs * sizeof(double *)); /* reply blocks storage */

    MPI_Status status;
    int req_count;
    for (int p = 0; p < all_process_expect[rank]; p++) {
        /* Wait until a request comes */
        MPI_Probe(MPI_ANY_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &req_count);
        rep_buf[p] = (double *) malloc_or_exit(req_count * sizeof(double));
        /* fill rep_buf[p] with requested x elements */
        MPI_Recv(reqs, req_count, MPI_INT, status.MPI_SOURCE, REQUEST_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i < req_count; i++) {
            if (reqs[i] < proc_info[rank].first_row || reqs[i] > proc_info[rank].last_row) {
                printf("Wrong index %d looking at process %d\n", reqs[i], p);
                return NULL;
            }
            rep_buf[p][i] = buf_x[reqs[i] - proc_info[rank].first_row];
        }

        /* send the requested block */
        MPI_Isend(rep_buf[p], req_count, MPI_DOUBLE, status.MPI_SOURCE, REPLY_TAG, MPI_COMM_WORLD, &send_reqs[0]);
//        printf("[%d] Replying requests from process %2d \t[%5d]\n", rank, status.MPI_SOURCE, req_count);
    }
//    printf("[%d] Replied to all requests! [%4d]\n", rank, to_send[rank]);

    /* Local elements multiplication */
    for (int k = 0; k < proc_info[rank].NZ; k++) {
        if (in_diagonal(buf_j_idx[k], proc_info[rank].first_row, proc_info[rank].last_row)) {
            y[buf_i_idx[k] - proc_info[rank].first_row] +=
                    buf_values[k] * buf_x[buf_j_idx[k] - proc_info[rank].first_row];
        }
    }



    /* wait for all blocks to arrive */
    int p;
//    printf("[%d] Waiting for %d requests\n", rank, req_made);
    double *vecFromRemotePros = (double *) calloc_or_exit(proc_info[rank].N, sizeof(double));
    for (int q = 0; q < req_made; q++) {
        MPI_Waitany(nprocs, recv_reqs, &p, MPI_STATUS_IGNORE);
        assert(p != MPI_UNDEFINED);

        /* fill x array with new elements */
        for (int i = 0; i < to_send[p]; i++) {
            vecFromRemotePros[send_buf[p][i]] = recv_buf[p][i];
        }
    }

    /* Global elements multiplication */
    for (int k = 0; k < proc_info[rank].NZ; k++) {
        if (!in_diagonal(buf_j_idx[k], proc_info[rank].first_row, proc_info[rank].last_row)) {
            y[buf_i_idx[k] - proc_info[rank].first_row] += buf_values[k] * vecFromRemotePros[buf_j_idx[k]];
        }
    }

    /* gather y elements from processes and save it to res */
    debug("[%d] Gathering results...\n", rank);
    MPI_Gatherv(y, proc_info[rank].M, MPI_DOUBLE, res, row_count,
                row_offset, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    /* return final result */
    return res;
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

int main(int argc, char *argv[]) {
    char *in_file,
            *out_file = NULL;

    double t, comp_time, partition_time;
    int nprocs,     /* number of tasks/processes */
            rank;       /* id of task/process */

    /***** MPI MASTER (root) process only ******/
//    proc_info_t *all_proc_info;

    int *buf_i_idx,     /* row index for all matrix elements */
            *buf_j_idx;     /* column index for all matrix elements */
    double *buf_values, /* value for all matrix elements */
            *buf_x,      /* value for all x vector elements */
            *res;        /* final result -> Ax */

    /*******************************************/

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    create_mpi_datatypes(&proc_info_type);

    /* master thread reads matrix */
    if (argc < 2 || argc > 3) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        in_file = argv[1];
        if (argc == 3)
            out_file = argv[2];
    }

    /* initialize proc_info array */
    proc_info = (proc_info_t *) malloc_or_exit(nprocs * sizeof(proc_info_t));

    if (rank_wise_read_matrix(in_file, &buf_i_idx, &buf_j_idx, &buf_values,
                              &proc_info[rank].M, &proc_info[rank].N, &proc_info[rank].NZ,
                              &proc_info[rank].first_row, &proc_info[rank].last_row, rank) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }

//    printf("[%d] Read matrix from '%s'!\n", rank, in_file);
//    printf("[%d] Matrix properties: M=%d, N = %d, NZ = %d, first_row=%d, last_row=%d\n\n", rank, proc_info[rank].M, proc_info[rank].N, proc_info[rank].NZ, proc_info[rank].first_row, proc_info[rank].last_row);

    for (int j = 0; j < proc_info[rank].NZ; ++j) {
//        printf("rank=%d, i=%d, j=%d, values=%lf\n", rank, buf_i_idx[j], buf_j_idx[j], buf_values[j]);
    }

    buf_x = (double *) malloc_or_exit(proc_info[rank].N * sizeof(double));
    res = (double *) malloc_or_exit(proc_info[rank].N * sizeof(double));
    for (int i = 0; i < proc_info[rank].N; i++) {
        buf_x[i] = 1;
    }
    if (rank == MASTER) {


        t = MPI_Wtime();

        partition_time = (MPI_Wtime() - t) * 1000.0;

        debug("[%d] Partition time: %10.3lf ms\n\n", rank, partition_time);
        debug("[%d] Starting algorithm...\n", rank);
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

    /* allocate buffers for requests sending */
    int **send_buf = (int **) malloc_or_exit(nprocs * sizeof(int *));
    for (int i = 0; i < nprocs; i++) {
        if (i != rank && proc_info[i].M > 0)
            send_buf[i] = (int *) malloc_or_exit(proc_info[i].M * sizeof(int));
    }

    int *to_send = (int *) calloc_or_exit(nprocs, sizeof(int));    /* # of req to each proc */
    int *map = (int *) calloc_or_exit(proc_info[rank].N, sizeof(int));

    /* build sending blocks to processors */
    int dest, col;
    for (int i = 0; i < proc_info[rank].NZ; i++) {
        col = buf_j_idx[i];
        /* check whether I need to send a request */
        if (in_diagonal(col, proc_info[rank].first_row, proc_info[rank].last_row) || map[col] > 0) {
            continue;
        }

        /* search which process has the element
         * NOTE: Due to small number or processes, serial search is faster
         */
        dest = -1;
        for (int p = 0; p < nprocs; p++) {
            if (in_diagonal(col, proc_info[p].first_row, proc_info[p].last_row)) {
                dest = p;
                break;
            }
        }
        assert(dest >= 0);
        /* insert new request */
        send_buf[dest][to_send[dest]++] = col;
        map[col] = 1;
    }

    /* Matrix-vector multiplication for each processes */
    res = mat_vec_mult_parallel(rank, nprocs, buf_i_idx, buf_j_idx, buf_values, buf_x, row_count, row_offset/*, send_buf,
                                to_send*/);
    if (rank == MASTER) {
        printf("Result Y= ");
        for (int i = 0; i < proc_info[MASTER].N; ++i) {
            printf("|%lf| ", res[i]);
        }
    }


    /*MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;*/
/*
    double stdev = 0, mean = 0, runs[TOTAL_RUNS];
    for (int r = 0; r < TOTAL_RUNS; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == MASTER) t = MPI_Wtime();
        res = mat_vec_mult_parallel(rank, nprocs, all_proc_info, buf_i_idx,
                                    buf_j_idx, buf_values, buf_x);
        //MPI_Barrier(MPI_COMM_WORLD);
        if (rank == MASTER){
            runs[r] = (MPI_Wtime() - t) * 1000.0;
            mean += runs[r];
        }
    }*/

    /* print execution stats */
    /*  if (rank == MASTER) {
          mean /= TOTAL_RUNS;
          for (int r = 0; r < TOTAL_RUNS; r++) {
              stdev += (runs[r] - mean) * (runs[r] - mean);
          }
          stdev = sqrt(stdev/(TOTAL_RUNS-1));

          printf("[%d] Computation time: %10.3lf [%4.3lf] ms\n\n", rank, mean, stdev);
          printf("[%d] Total execution time: %10.3lf ms\n", rank, mean + partition_time);
          debug("Finished!\n");
      }
  */
    /* write to output file */
    if (rank == MASTER) {

        /*FILE *resultCSV;
        FILE *checkFile;
        if((checkFile = fopen("MPISpMVResult.csv","r"))!=NULL)
        {
            // file exists
            fclose(checkFile);
            if ( !(resultCSV = fopen("MPISpMVResult.csv", "a")) ) {
                fprintf(stderr, "fopen: failed to open file MPISpMVResult.csv");
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            if ( !(resultCSV = fopen("MPISpMVResult.csv", "w")) ) {
                fprintf(stderr, "fopen: failed to open file MPISpMVResult.csv");
                exit(EXIT_FAILURE);
            }
            fprintf(resultCSV, "MatrixName,ComputationTime,Stdev,TotalRun,nProcess,PartitionType,TotalExecutionTime\n");
        }

        fprintf(resultCSV, "%s,%10.3lf,%4.3lf,%d,%d,%d,%10.3lf\n", in_file, mean, stdev, TOTAL_RUNS, nprocs, (policy==EQUAL_ROWS?0:1), (mean + partition_time));
        if ( fclose(resultCSV) != 0) {
            fprintf(stderr, "fopen: failed to open file MPISpMVResult");
            exit(EXIT_FAILURE);
        }
        if (out_file != NULL) {
            printf("Writing result to '%s'\n", out_file);

            *//* open file *//*
            FILE *f;
            if ( !(f = fopen(out_file, "w")) ) {
                fprintf(stderr, "fopen: failed to open file '%s'", out_file);
                exit(EXIT_FAILURE);
            }

            fprintf(f, "Vector:\n");
            for (int i = 0; i < all_proc_info[MASTER].N; i++) {
                fprintf(f, "%.8lf ", buf_x[i]);
            }
            fprintf(f, "\n Result:\n");
            *//* write result *//*
            for (int i = 0; i < all_proc_info[MASTER].N; i++) {
                fprintf(f, "%.8lf\n", res[i]);
            }

            *//* close file *//*
            if ( fclose(f) != 0) {
                fprintf(stderr, "fopen: failed to open file '%s'", out_file);
                exit(EXIT_FAILURE);
            }

            printf("Done!\n");
        }*/

        /* free the memory */
        /*free(buf_values);
        free(buf_i_idx);
        free(buf_j_idx);
        free(buf_x);
        free(res);*/
    }
    free(buf_values);
    free(buf_i_idx);
    free(buf_j_idx);
    free(buf_x);
    free(res);

    /* MPI: end */
    MPI_Finalize();

    return 0;
}


/* 
 * High-level wrapper for Matrix Market I/O library
 *
 * It is used to read sparse, real & square matrices from files
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mmio-wrapper.h"
#include "mmio.h"
#include "util.h"
#include <math.h>

/* Reads a matrix from a Matrix Market file, stored in COO format */
int read_matrix(const char *filename, int **i_idx, int **j_idx, double **values, int *N, int *NZ) {
    FILE *f;
    MM_typecode matcode;
    int errorcode, nrows, ncols, nz_elements;

    /* open the file */
    if ((f = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Cannot open '%s'\n", filename);
        return 1;
    }

    /* process first line */
    if ((errorcode = mm_read_banner(f, &matcode)) != 0) {
        fprintf(stderr, "Error while processing banner (file:'%s') (code=%d)\n",
                filename, errorcode);
        return 1;
    }

    /* matrix should be sparse and real */
    if (!mm_is_matrix(matcode) ||
        !mm_is_real(matcode) ||
        !mm_is_sparse(matcode)) {
        fprintf(stderr, "Not supported matrix type: %s\n", mm_typecode_to_str(matcode));
        return 1;
    }

    /* read info */
    if ((errorcode = mm_read_mtx_crd_size(f, &nrows, &ncols, &nz_elements)) != 0) {
        fprintf(stderr, "Error while processing array (file:'%s') (code:%d)\n",
                filename, errorcode);
        return 1;
    }

    /* matrix should be square */
    if (nrows != ncols) {
        fprintf(stderr, "Matrix is NOT square (rows=%d, cols=%d)\n", nrows, ncols);
        return 1;
    }

    *N = nrows;
    *NZ = nz_elements;

    /* reserve memory for vector */
    *i_idx = (int *) malloc(nz_elements * sizeof(int));
    *j_idx = (int *) malloc(nz_elements * sizeof(int));
    *values = (double *) malloc(nz_elements * sizeof(double));

    /* read actual matrix */
    for (int i = 0; i < *NZ; i++) {
        fscanf(f, "%d %d %lf", &(*i_idx)[i], &(*j_idx)[i], &(*values)[i]);
        (*i_idx)[i]--;
        (*j_idx)[i]--;
    }

    /* close the file */
    if (fclose(f) != 0) {
        fprintf(stderr, "Cannot close file (fil:'%s')\n", filename);
    }

    return 0;
}

int rank_wise_read_matrix(const char *filename, int **i_idx, int **j_idx, double **values, int *M, int *N, int *NZ,
                          int *first_row, int *last_row, int rank) {
    FILE *f;
    MM_typecode matcode;
    int errorcode, nrows, ncols, nz_elements, start_row = 0, end_row = 0;

    /* open the file */
    char rank_wise_filename[MM_MAX_LINE_LENGTH];
    char *_ptr = strtok(filename, ".");
    sprintf(rank_wise_filename, "%s_%d.%s", _ptr, rank, strtok(NULL, "-"));
    if ((f = fopen(rank_wise_filename, "r")) == NULL) {
        fprintf(stderr, "Cannot open '%s'\n", rank_wise_filename);
        return 1;
    }

    /* process first line */
    if ((errorcode = mm_read_banner(f, &matcode)) != 0) {
        fprintf(stderr, "Error while processing banner (file:'%s') (code=%d)\n",
                filename, errorcode);
        return 1;
    }

    /* matrix should be sparse and real */
    if (!mm_is_matrix(matcode) ||
        !mm_is_real(matcode) ||
        !mm_is_sparse(matcode)) {
        fprintf(stderr, "Not supported matrix type: %s\n", mm_typecode_to_str(matcode));
        return 1;
    }

    /* read info */
    if ((errorcode = mm_read_mtx_crd_size(f, &nrows, &ncols, &nz_elements)) != 0) {
        fprintf(stderr, "Error while processing array (file:'%s') (code:%d)\n",
                filename, errorcode);
        return 1;
    }

    /* matrix should be square */
    /*if (nrows != ncols) {
        fprintf(stderr, "Matrix is NOT square (rows=%d, cols=%d)\n", nrows, ncols);
        return 1;
    }*/

    *M = nrows;
    *N = ncols;
    *NZ = nz_elements;
    start_row = ncols;

    /* reserve memory for vector */
    *i_idx = (int *) malloc(nz_elements * sizeof(int));
    *j_idx = (int *) malloc(nz_elements * sizeof(int));
    *values = (double *) malloc(nz_elements * sizeof(double));

    /* read actual matrix */
    for (int i = 0; i < *NZ; i++) {
        fscanf(f, "%d %d %lf", &(*i_idx)[i], &(*j_idx)[i], &(*values)[i]);
        (*i_idx)[i]--;
        (*j_idx)[i]--;
        if (start_row > (*i_idx)[i])
            start_row = (*i_idx)[i];
        if (end_row < (*i_idx)[i])
            end_row = (*i_idx)[i];
    }

    (*first_row) = start_row;
    (*last_row) = end_row;
    /* close the file */
    if (fclose(f) != 0) {
        fprintf(stderr, "Cannot close file (fil:'%s')\n", filename);
    }

    return 0;
}

int test_csr_read_2D_partitioned_mat(const char *filename, int **row_ptr, int **col_ptr, double **val_ptr,
                                int sqrRank, int rank) {
    FILE *f;
    MM_typecode matcode;
    int errorcode, nrows, ncols, nz_elements;

    /* open the file */
    char rank_wise_filename[MM_MAX_LINE_LENGTH];
    char *_ptr = strtok(filename, ".");
    sprintf(rank_wise_filename, "%s_%d.%s", _ptr, rank, strtok(NULL, "-"));
    if ((f = fopen(rank_wise_filename, "r")) == NULL) {
        printf("Cannot open '%s'\n", rank_wise_filename);
        return 1;
    }

    /* process first line */
    if ((errorcode = mm_read_banner(f, &matcode)) != 0) {
        printf("Error while processing banner (file:'%s') (code=%d)\n", filename, errorcode);
        return 1;
    }

    /* matrix should be sparse and real */
    if (!mm_is_matrix(matcode) ||
        !mm_is_real(matcode) ||
        !mm_is_sparse(matcode)) {
        printf("Not supported matrix type: %s\n", mm_typecode_to_str(matcode));
        return 1;
    }

    /* read info */
    if ((errorcode = mm_read_mtx_crd_size(f, &nrows, &ncols, &nz_elements)) != 0) {
        printf("Error while processing array (file:'%s') (code:%d)\n", filename, errorcode);
        return 1;
    }

    int startRow = ceil(((double) ncols / sqrRank)) * (rank / sqrRank);
    if (nrows <= 0) {
        printf("[%d] issue with rows=%d\n", rank, nrows);
    }
    if (ncols <= 0) {
        printf("[%d] issue with columns=%d\n", rank, ncols);
    }
    if (nz_elements <= 0) {
        printf("[%d] issue with nz_elements=%d\n", rank, nz_elements);
    }
    /// Initialize CSR row, col and value pointer.
    (*row_ptr) = (int *) calloc_or_exit((nrows + 1), sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(nz_elements * sizeof(int));
    (*val_ptr) = (double *) malloc_or_exit(nz_elements * sizeof(double));

    (*row_ptr)[0] = 0;
    int *i_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    int *j_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    double *values = (double *) malloc_or_exit(nz_elements * sizeof(double));
    /* read actual matrix */
    for (int i = 0; i < nz_elements; i++) {
        fscanf(f, "%d %d %lf", &(i_idx[i]), &(j_idx[i]), &(values[i]));
        i_idx[i]--;
        j_idx[i]--;
    }
    for (int i = 0; i < nz_elements; i++) {
        if ((i_idx[i] - startRow) >= nrows || (i_idx[i] - startRow) < 0) {
            printf("[%d] Index out of bound for row=%d, start row=%d\n", rank, i_idx[i], startRow);
            return 1;
        }
        (*row_ptr)[i_idx[i] - startRow]++;
    }
    for (int i = 0, cumsum = 0; i < nrows; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = cumsum;
        cumsum += temp;
    }
    (*row_ptr)[nrows] = nz_elements;

    for (int n = 0; n < nz_elements; n++) {
        int row = i_idx[n] - startRow;
        if (row < 0 || row >= nrows) {
            printf("[%d] out of bound for row=%d, start row=%d\n", rank, row, startRow);
            return 1;
        }
        int dest = (*row_ptr)[row];
        (*col_ptr)[dest] = j_idx[n];
        (*val_ptr)[dest] = values[n];

        (*row_ptr)[row]++;
    }

    for (int i = 0, last = 0; i <= nrows; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = last;
        last = temp;
    }
    /* close the file */
    if (fclose(f) != 0) {
        fprintf(stderr, "Cannot close file (fil:'%s')\n", filename);
    }

    return 0;
}

int csr_read_2D_partitioned_mat(const char *filename, int **row_ptr, int **col_ptr, double **val_ptr,
                                proc_info_t **ranks_info, int sqrRank, int rank) {
    FILE *f;
    MM_typecode matcode;
    int errorcode, nrows, ncols, nz_elements;

    /* open the file */
    char rank_wise_filename[MM_MAX_LINE_LENGTH];
    char *_ptr = strtok(filename, ".");
    sprintf(rank_wise_filename, "%s_%d.%s", _ptr, rank, strtok(NULL, "-"));
    if ((f = fopen(rank_wise_filename, "r")) == NULL) {
        printf("Cannot open '%s'\n", rank_wise_filename);
        return 1;
    }

    /* process first line */
    if ((errorcode = mm_read_banner(f, &matcode)) != 0) {
        printf("Error while processing banner (file:'%s') (code=%d)\n", filename, errorcode);
        return 1;
    }

    /* matrix should be sparse and real */
    if (!mm_is_matrix(matcode) ||
        !mm_is_real(matcode) ||
        !mm_is_sparse(matcode)) {
        printf("Not supported matrix type: %s\n", mm_typecode_to_str(matcode));
        return 1;
    }

    /* read info */
    if ((errorcode = mm_read_mtx_crd_size(f, &nrows, &ncols, &nz_elements)) != 0) {
        printf("Error while processing array (file:'%s') (code:%d)\n", filename, errorcode);
        return 1;
    }

    int startRow = ceil(((double) ncols / sqrRank)) * (rank / sqrRank);
    (*ranks_info)[rank].M = ceil(((double)ncols)/sqrRank);
    (*ranks_info)[rank].N = ncols;
    (*ranks_info)[rank].NZ = nz_elements;
    (*ranks_info)[rank].first_row = startRow;
    (*ranks_info)[rank].last_row = startRow + (*ranks_info)[rank].M - 1;
    /// Initialize CSR row, col and value pointer.
    (*row_ptr) = (int *) calloc_or_exit(((*ranks_info)[rank].M + 1), sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(nz_elements * sizeof(int));
    (*val_ptr) = (double *) malloc_or_exit(nz_elements * sizeof(double));

    (*row_ptr)[0] = 0;
    int *i_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    int *j_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    double *values = (double *) malloc_or_exit(nz_elements * sizeof(double));
    /* read actual matrix */
    for (int i = 0; i < nz_elements; i++) {
        fscanf(f, "%d %d %lf", &(i_idx[i]), &(j_idx[i]), &(values[i]));
        i_idx[i]--;
        j_idx[i]--;
    }
    for (int i = 0; i < nz_elements; i++) {
        if ((i_idx[i] - startRow) >= (*ranks_info)[rank].M || (i_idx[i] - startRow) < 0) {
            printf("[%d] Index out of bound for row=%d, start row=%d\n", rank, i_idx[i], startRow);
            return 1;
        }
        (*row_ptr)[i_idx[i] - startRow]++;
    }
    for (int i = 0, cumsum = 0; i < (*ranks_info)[rank].M; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = cumsum;
        cumsum += temp;
    }
    (*row_ptr)[(*ranks_info)[rank].M] = nz_elements;

    for (int n = 0; n < nz_elements; n++) {
        int row = i_idx[n] - startRow;
        if (row < 0 || row >= (*ranks_info)[rank].M) {
            printf("[%d] out of bound for row=%d, start row=%d\n", rank, row, startRow);
            return 1;
        }
        int dest = (*row_ptr)[row];
        (*col_ptr)[dest] = j_idx[n];
        (*val_ptr)[dest] = values[n];

        (*row_ptr)[row]++;
    }

    for (int i = 0, last = 0; i <= (*ranks_info)[rank].M; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = last;
        last = temp;
    }
    /* close the file */
    if (fclose(f) != 0) {
        fprintf(stderr, "Cannot close file (fil:'%s')\n", filename);
    }

    return 0;
}

int rank_wise_read_matrix_csr(const char *filename, int **row_ptr, int **col_ptr, double **val_ptr,
                              proc_info_t **ranks_info, int rank, int *offDiagonalElements) {
    FILE *f;
    MM_typecode matcode;
    int errorcode, nrows, ncols, nz_elements, start_row = 0, end_row = 0;

    /* open the file */
    char rank_wise_filename[MM_MAX_LINE_LENGTH];
    char *_ptr = strtok(filename, ".");
    sprintf(rank_wise_filename, "%s_%d.%s", _ptr, rank, strtok(NULL, "-"));
    if ((f = fopen(rank_wise_filename, "r")) == NULL) {
        fprintf(stderr, "Cannot open '%s'\n", rank_wise_filename);
        return 1;
    }

    /* process first line */
    if ((errorcode = mm_read_banner(f, &matcode)) != 0) {
        fprintf(stderr, "Error while processing banner (file:'%s') (code=%d)\n",
                filename, errorcode);
        return 1;
    }

    /* matrix should be sparse and real */
    if (!mm_is_matrix(matcode) ||
        !mm_is_real(matcode) ||
        !mm_is_sparse(matcode)) {
        fprintf(stderr, "Not supported matrix type: %s\n", mm_typecode_to_str(matcode));
        return 1;
    }

    /* read info */
    if ((errorcode = mm_read_mtx_crd_size(f, &nrows, &ncols, &nz_elements)) != 0) {
        fprintf(stderr, "Error while processing array (file:'%s') (code:%d)\n",
                filename, errorcode);
        return 1;
    }
    (*ranks_info)[rank].M = nrows;
    (*ranks_info)[rank].N = ncols;
    (*ranks_info)[rank].NZ = nz_elements;

    /// Initialize CSR row, col and value pointer.
    (*row_ptr) = (int *) calloc_or_exit((nrows + 1), sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(nz_elements * sizeof(int));
    (*val_ptr) = (double *) malloc_or_exit(nz_elements * sizeof(double));

    start_row = ncols;
    (*row_ptr)[0] = 0;
    int *i_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    int *j_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    double *values = (double *) malloc_or_exit(nz_elements * sizeof(double));
    /* read actual matrix */
    for (int i = 0; i < nz_elements; i++) {
        fscanf(f, "%d %d %lf", &(i_idx[i]), &(j_idx[i]), &(values[i]));
        i_idx[i]--;
        j_idx[i]--;
        if (start_row > i_idx[i])
            start_row = i_idx[i];
        if (end_row < i_idx[i])
            end_row = i_idx[i];
    }
    (*ranks_info)[rank].first_row = start_row;
    (*ranks_info)[rank].last_row = end_row;
    for (int i = 0; i < nz_elements; i++) {
        if ((i_idx[i] - start_row) >= nrows || (i_idx[i] - start_row) < 0) {
            printf("[%d] Index out of bound for row=%d, start row=%d, end row=%d\n", rank, i_idx[i], start_row,
                   end_row);
        }
        (*row_ptr)[i_idx[i] - start_row]++;
    }

    for (int i = 0, cumsum = 0; i < nrows; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = cumsum;
        cumsum += temp;
    }
    (*row_ptr)[nrows] = nz_elements;

    for (int n = 0; n < nz_elements; n++) {
        int row = i_idx[n] - start_row;
        if (row < 0 || row >= nrows) {
            printf("[%d] out of bound for row=%d, start row=%d, end row=%d\n", rank, row, start_row, end_row);
        }
        int dest = (*row_ptr)[row];
        if (!in_diagonal(j_idx[n], start_row, end_row)) {
            (*offDiagonalElements)++;
        }
        (*col_ptr)[dest] = j_idx[n];
        (*val_ptr)[dest] = values[n];

        (*row_ptr)[row]++;
    }

    for (int i = 0, last = 0; i <= nrows; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = last;
        last = temp;
    }
    /* close the file */
    if (fclose(f) != 0) {
        fprintf(stderr, "Cannot close file (fil:'%s')\n", filename);
    }

    return 0;
}

int write_matrix(const char *filename, const int *i_idx, const int *j_idx, const double *values, int N, int NZ) {
    FILE *f;
    MM_typecode matcode;

    /* open the file */
    if ((f = fopen(filename, "w")) == NULL) {
        fprintf(stderr, "Cannot open '%s'\n", filename);
        return 1;
    }

    /* init and set proper flags for matrix */
    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_real(&matcode);
    mm_set_sparse(&matcode);

    /* write banner and matrix size info */
    mm_write_banner(f, matcode);
    mm_write_mtx_crd_size(f, N, N, NZ);

    /* write matrix elements */
    for (int i = 0; i < NZ; i++) {
        fprintf(f, "%d %d %.9lf\n", i_idx[i] + 1, j_idx[i] + 1, values[i]);
    }

    /* close the file */
    if (fclose(f) != 0) {
        fprintf(stderr, "Cannot close file (fil:'%s')\n", filename);
    }

    return 0;
}


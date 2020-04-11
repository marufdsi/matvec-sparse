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
    char file_n[100];
    strcpy(file_n, filename);
    char *_ptr = strtok(file_n, ".");
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
    char file_n[100];
    strcpy(file_n, filename);
    char *_ptr = strtok(file_n, ".");
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
    nrows = 123761;
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

int csr_read_2D_partitioned_mat(const char *filename, int **row_ptr, int **col_ptr, f_type **val_ptr,
                                proc_info_t **ranks_info, int sqrRank, int rank) {
    FILE *f;
    MM_typecode matcode;
    int errorcode, nrows, ncols, nz_elements;

    /* open the file */
    char rank_wise_filename[MM_MAX_LINE_LENGTH];
    char file_n[100];
    strcpy(file_n, filename);
//    char *_ptr = strtok(file_n, ".");
//    sprintf(rank_wise_filename, "%s_%d.%s", _ptr, rank, strtok(NULL, "-"));
    printf("file name: %s\n", file_n);
    char *_ptr = strtok(file_n, ".mtx");
    sprintf(rank_wise_filename, "%s_%d.mtx", _ptr, rank);
    printf("rank-wise file name: %s\n", rank_wise_filename);
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
    (*val_ptr) = (f_type *) malloc_or_exit(nz_elements * sizeof(f_type));

    (*row_ptr)[0] = 0;
    int *i_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    int *j_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    f_type *values = (f_type *) malloc_or_exit(nz_elements * sizeof(f_type));
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
    char file_n[100];
    strcpy(file_n, filename);
    char *_ptr = strtok(file_n, ".");
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

int read_coo_matrix_to_csr(const char *filename, int **row_ptr, int **col_ptr, ValueType **val_ptr, int *mat_row, int *_nnz) {
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
    (*mat_row) = nrows;
    (*_nnz) = nz_elements;

    /// Initialize CSR row, col and value pointer.
    (*row_ptr) = (int *) calloc_or_exit((nrows + 1), sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(nz_elements * sizeof(int));
    (*val_ptr) = (ValueType *) malloc_or_exit(nz_elements * sizeof(ValueType));

    (*row_ptr)[0] = 0;
    int *i_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    int *j_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    ValueType *values = (ValueType *) malloc_or_exit(nz_elements * sizeof(ValueType));
    /* read actual matrix */
    for (int i = 0; i < nz_elements; i++) {
        fscanf(f, "%d %d %lf", &(i_idx[i]), &(j_idx[i]), &(values[i]));
        i_idx[i]--;
        j_idx[i]--;
    }
    for (int i = 0; i < nz_elements; i++) {
        if ((i_idx[i]) >= nrows || (i_idx[i]) < 0) {
            printf("Index out of bound for row=%d\n", i_idx[i]);
        }
        (*row_ptr)[i_idx[i]]++;
    }

    for (int i = 0, cumsum = 0; i < nrows; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = cumsum;
        cumsum += temp;
    }
    (*row_ptr)[nrows] = nz_elements;

    for (int n = 0; n < nz_elements; n++) {
        int row = i_idx[n];
        if (row < 0 || row >= nrows) {
            printf("out of bound for row=%d\n", row);
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

/*int csr_read_Matrix(const char *filename, int **row_ptr, int **col_ptr, ValueType **val_ptr, int *mat_row, int *mat_col, int *_nnz, int *max_deg){
    // load matrix
    FILE *f;
    MM_typecode matcode;
    int isInteger = 0, isReal = 0, isPattern = 0, isSymmetric = 0, ret_code, m, n, nnzA, nnzA_mtx_report;
    if ((f = fopen(filename, "r")) == NULL)
        return -1;

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        return -2;
    }

    if ( mm_is_complex( matcode ) )
    {
        printf("Sorry, data type 'COMPLEX' is not supported.\n");
        return -3;
    }

    if ( mm_is_pattern( matcode ) )  { isPattern = 1; *//*cout << "type = Pattern" << endl;*//* }
    if ( mm_is_real ( matcode) )     { isReal = 1; *//*cout << "type = real" << endl;*//* }
    if ( mm_is_integer ( matcode ) ) { isInteger = 1; *//*cout << "type = integer" << endl;*//* }

    *//* find out size of sparse matrix .... *//*
    ret_code = mm_read_mtx_crd_size(f, &m, &n, &nnzA_mtx_report);
    if (ret_code != 0)
        return -4;

    (*mat_row) = m;
    (*mat_col) = n;
    (*_nnz) = nnzA_mtx_report;
    if ( mm_is_symmetric( matcode ) || mm_is_hermitian( matcode ) )
    {
        isSymmetric = 1;
        //cout << "symmetric = true" << endl;
    }
    else
    {
        //cout << "symmetric = false" << endl;
    }

    int *csrRowPtrA_counter = (int *)malloc((m+1) * sizeof(int));
    memset(csrRowPtrA_counter, 0, (m+1) * sizeof(int));

    int *csrRowIdxA_tmp = (int *)malloc(nnzA_mtx_report * sizeof(int));
    int *csrColIdxA_tmp = (int *)malloc(nnzA_mtx_report * sizeof(int));
    ValueType *csrValA_tmp    = (ValueType *)malloc(nnzA_mtx_report * sizeof(ValueType));

    *//* NOTE: when reading in doubles, ANSI C requires the use of the "l"  *//*
    *//*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur *//*
    *//*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            *//*

    for (int i = 0; i < nnzA_mtx_report; i++)
    {
        int idxi, idxj;
        double fval;
        int ival;

        if (isReal)
            fscanf(f, "%d %d %lg\n", &idxi, &idxj, &fval);
        else if (isInteger)
        {
            int count = fscanf(f, "%d %d %d\n", &idxi, &idxj, &ival);
            fval = ival;
        }
        else if (isPattern)
        {
            int count = fscanf(f, "%d %d\n", &idxi, &idxj);
            fval = 1.0;
        }

        // adjust from 1-based to 0-based
        idxi--;
        idxj--;

        csrRowPtrA_counter[idxi]++;
        csrRowIdxA_tmp[i] = idxi;
        csrColIdxA_tmp[i] = idxj;
        csrValA_tmp[i] = fval;
        if((*max_deg) < csrRowPtrA_counter[idxi]){
            (*max_deg) = csrRowPtrA_counter[idxi];
        }
    }

    if (f != stdin)
        fclose(f);

    if (isSymmetric)
    {
        for (int i = 0; i < nnzA_mtx_report; i++)
        {
            if (csrRowIdxA_tmp[i] != csrColIdxA_tmp[i])
                csrRowPtrA_counter[csrColIdxA_tmp[i]]++;
        }
    }

    // exclusive scan for csrRowPtrA_counter
    int old_val, new_val;

    old_val = csrRowPtrA_counter[0];
    csrRowPtrA_counter[0] = 0;
    for (int i = 1; i <= m; i++)
    {
        new_val = csrRowPtrA_counter[i];
        csrRowPtrA_counter[i] = old_val + csrRowPtrA_counter[i-1];
        old_val = new_val;
    }

    nnzA = csrRowPtrA_counter[m];
    (*row_ptr) = (int *)malloc((m+1) * sizeof(int));
    memcpy((*row_ptr), csrRowPtrA_counter, (m+1) * sizeof(int));
    memset(csrRowPtrA_counter, 0, (m+1) * sizeof(int));

    (*col_ptr) = (int *)malloc(nnzA * sizeof(int));
    (*val_ptr)    = (ValueType *)malloc(nnzA * sizeof(ValueType));

    double gb = (double)((m + 1 + nnzA) * sizeof(int) + (2 * nnzA + m) * sizeof(ValueType));
    double gflop = (double)(2 * nnzA);

    if (isSymmetric)
    {
        for (int i = 0; i < nnzA_mtx_report; i++)
        {
            if (csrRowIdxA_tmp[i] != csrColIdxA_tmp[i])
            {
                int offset = (*row_ptr)[csrRowIdxA_tmp[i]] + csrRowPtrA_counter[csrRowIdxA_tmp[i]];
                (*col_ptr)[offset] = csrColIdxA_tmp[i];
                (*val_ptr)[offset] = csrValA_tmp[i];
                csrRowPtrA_counter[csrRowIdxA_tmp[i]]++;

                offset = (*row_ptr)[csrColIdxA_tmp[i]] + csrRowPtrA_counter[csrColIdxA_tmp[i]];
                (*col_ptr)[offset] = csrRowIdxA_tmp[i];
                (*val_ptr)[offset] = csrValA_tmp[i];
                csrRowPtrA_counter[csrColIdxA_tmp[i]]++;
            }
            else
            {
                int offset = (*row_ptr)[csrRowIdxA_tmp[i]] + csrRowPtrA_counter[csrRowIdxA_tmp[i]];
                (*col_ptr)[offset] = csrColIdxA_tmp[i];
                (*val_ptr)[offset] = csrValA_tmp[i];
                csrRowPtrA_counter[csrRowIdxA_tmp[i]]++;
            }
        }
    }
    else
    {
        for (int i = 0; i < nnzA_mtx_report; i++)
        {
            int offset = (*row_ptr)[csrRowIdxA_tmp[i]] + csrRowPtrA_counter[csrRowIdxA_tmp[i]];
            (*col_ptr)[offset] = csrColIdxA_tmp[i];
            (*val_ptr)[offset] = csrValA_tmp[i];
            csrRowPtrA_counter[csrRowIdxA_tmp[i]]++;
        }
    }

    // free tmp space
    free(csrColIdxA_tmp);
    free(csrValA_tmp);
    free(csrRowIdxA_tmp);
    free(csrRowPtrA_counter);
    return 0;
}*/

/*int read_coo_matrix_to_csr_with_max_deg(const char *filename, int **row_ptr, int **col_ptr, ValueType **val_ptr, int *mat_row, int *_nnz, int *max_deg) {
    FILE *f;
    MM_typecode matcode;
    int errorcode, nrows, ncols, nz_elements;

    *//* open the file *//*
    if ((f = fopen(filename, "r")) == NULL) {
        fprintf(stderr, "Cannot open '%s'\n", filename);
        return 1;
    }

    *//* process first line *//*
    if ((errorcode = mm_read_banner(f, &matcode)) != 0) {
        fprintf(stderr, "Error while processing banner (file:'%s') (code=%d)\n",
                filename, errorcode);
        return 1;
    }

    *//* matrix should be sparse and real *//*
    if (!mm_is_matrix(matcode) ||
        !mm_is_real(matcode) ||
        !mm_is_sparse(matcode)) {
        fprintf(stderr, "Not supported matrix type: %s\n", mm_typecode_to_str(matcode));
        return 1;
    }

    *//* read info *//*
    if ((errorcode = mm_read_mtx_crd_size(f, &nrows, &ncols, &nz_elements)) != 0) {
        fprintf(stderr, "Error while processing array (file:'%s') (code:%d)\n",
                filename, errorcode);
        return 1;
    }
    (*mat_row) = nrows;
    (*_nnz) = nz_elements;
    (*max_deg) = 0;
    /// Initialize CSR row, col and value pointer.
    (*row_ptr) = (int *) calloc_or_exit((nrows + 1), sizeof(int));
    (*col_ptr) = (int *) malloc_or_exit(nz_elements * sizeof(int));
    (*val_ptr) = (ValueType *) malloc_or_exit(nz_elements * sizeof(ValueType));

    (*row_ptr)[0] = 0;
    int *i_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    int *j_idx = (int *) malloc_or_exit(nz_elements * sizeof(int));
    ValueType *values = (ValueType *) malloc_or_exit(nz_elements * sizeof(ValueType));
    *//* read actual matrix *//*
    for (int i = 0; i < nz_elements; i++) {
        fscanf(f, "%d %d %lf", &(i_idx[i]), &(j_idx[i]), &(values[i]));
        i_idx[i]--;
        j_idx[i]--;
        (*row_ptr)[i_idx[i]]++;
        if((*max_deg) < (*row_ptr)[i_idx[i]])
            (*max_deg) = (*row_ptr)[i_idx[i]];
    }
    *//*for (int i = 0; i < nz_elements; i++) {
        if ((i_idx[i]) >= nrows || (i_idx[i]) < 0) {
            printf("Index out of bound for row=%d\n", i_idx[i]);
        }
        (*row_ptr)[i_idx[i]]++;
        if(max_deg < (*row_ptr)[i_idx[i]])
            max_deg = (*row_ptr)[i_idx[i]];
    }*//*

    for (int i = 0, cumsum = 0; i < nrows; i++) {
        int temp = (*row_ptr)[i];
        (*row_ptr)[i] = cumsum;
        cumsum += temp;
    }
    (*row_ptr)[nrows] = nz_elements;

    for (int n = 0; n < nz_elements; n++) {
        int row = i_idx[n];
        if (row < 0 || row >= nrows) {
            printf("out of bound for row=%d\n", row);
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
    *//* close the file *//*
    if (fclose(f) != 0) {
        fprintf(stderr, "Cannot close file (fil:'%s')\n", filename);
    }

    return 0;
}*/

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


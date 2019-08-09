/* 
 * High-level wrapper for Matrix Market I/O library
 *
 * It is used to read sparse, real & square matrices from files
 */

#ifndef MM_IO_WRAPPER_H
#define MM_IO_WRAPPER_H

#include "util.h"

int read_matrix (const char * filename, int **i_idx, int **j_idx, double **values, int *N, int *NZ);
int rank_wise_read_matrix (const char * filename, int **i_idx, int **j_idx, double **values, int *M, int *N, int *NZ, int *first_row, int *last_row, int rank);
int rank_wise_read_matrix_csr (const char * filename, int **row_ptr, int **col_ptr, double **val_ptr, proc_info_t **ranks_info, int rank, int *offDiagonalElements);
int write_matrix (const char *filename, const int *i_idx, const int *j_idx, const double *values, int N, int NZ);

#endif


#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>

/* define useful debugging print macro */
#ifdef DEBUG
    #define debug(...) fprintf(stderr, __VA_ARGS__);
#else
    #define debug(...) do ; while(0)
#endif

typedef struct {
    int M,              /* Number of row of the processor */
        N,              /* dim of matrix */
        NZ,             /* number of non-zero elements */
        nz_count,       /* number of matrix elements for each process */
        nz_start_idx,   /* first matrix element for each process */
        row_count,      /* number of rows to process for each process */
        row_start_idx,  /* first row processed for each process */
        first_row,      /* Start row of the process */
        last_row;        /* End row of the process */
} proc_info_t;

struct Key {
    int col;
};
struct Value {
    double val;
};
typedef struct {
    struct Key key;
    struct Value value;
} Map;

double getVal(Map *map, int col, int limit);

void random_vec(double *v, int N, int limit);
int random_mat (int *buf_i, int *buf_j, double *buf_val, int start_row, int mat_row, int nzPerRow);
int csr_random_mat (int rank, proc_info_t *procs_info, int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int mat_col, int nzPerRow, int dist);
int csr_random_diagonal_mat (int rank, int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int nzPerRow);

void * malloc_or_exit(size_t size);
void * calloc_or_exit(size_t nmemb, size_t size);

char in_range(int num, int start, int count);
char in_diagonal(int num, int start, int end);

#endif

#include <time.h>

#include "util.h"

double getVal(Map *map, int col) {
    while (map != NULL) {
        if ((*map).key.col == col) {
            return (*map).value.val;
        }
        map++;
    }
    printf("Error!!! column=%d not found\n", col);
    return 0;
}

/*
 * Generate a random real vector of size N with
 * elements beloning to [0, limit)
 */
void random_vec (double *v, int N, int limit)
{
    /* fill v with random doubles */
    limit--;
    srand( 410 );
    for (int i = 0; i < N; i++) {
        v[i] = ((double)rand()) / (((double)RAND_MAX) / limit) + \
               ((double)rand()) / ((double)RAND_MAX);
    }
}

int random_mat (int *buf_i, int *buf_j, double *buf_val, int start_row, int mat_row, int nzPerRow)
{
    int start_idx = 0;
    for (int r = 0; r < mat_row; ++r) {
        srand(time(0));
        int *checkRepeat = (int *) calloc_or_exit(mat_row, sizeof(int));
        for (int i = 0; i < nzPerRow; i++) {
//            int conflict = 0;
            int rand_idx;
            do{
//                conflict++;
                rand_idx = rand() % mat_row;
            }while(checkRepeat[rand_idx]>0);
            /*if(conflict>0){
                printf("Conflicts Happened = %d\n", conflict);
            }*/
            checkRepeat[rand_idx] = 1;
            int idx = start_row + rand_idx;
            buf_i[start_idx] = start_row + r;
            buf_j[start_idx] = idx;
            /// Fill by any random double value
            buf_val[start_idx] = (double)(1 + (idx %10));
            start_idx++;
        }
    }

    return 1;
}

int csr_random_mat (int rank, int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int mat_col, int nzPerRow)
{
    int start_idx = 0;
    int row_elements = 0;
    row_ptr[0] = row_elements;
    srand(time(NULL));
    for (int r = 0; r < mat_row; ++r) {
        row_elements += nzPerRow;
        Map *map = (Map *) malloc_or_exit(nzPerRow * sizeof(Map));
        int off_diagonal = 0;
        int range = mat_col;
        int range_start = 0;
        for (int i = 0; i < nzPerRow; i++) {
            /*if((nzPerRow*2)/100 <= off_diagonal){
                range = mat_row;
                range_start = rank * mat_row;
            }*/
            int rand_idx;
            /// escape same random column
            do{
                rand_idx = rand() % range;
            }while(getVal(map, rand_idx)>0);
            if (rand_idx>= mat_row)
                off_diagonal++;
            map[i].key.col = rand_idx;
            map[i].value.val = 1.0;
            col_ptr[start_idx] = range_start + rand_idx;
            /// Fill by any random double value
            val_ptr[start_idx] = (double)(1 + (rand_idx %10));
            start_idx++;
        }
        row_ptr[r+1] = row_elements;
    }

    return 1;
}

int csr_random_diagonal_mat (int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int nzPerRow)
{
    int start_idx = 0;
    int row_elements = 0;
    row_ptr[0] = row_elements;
    for (int r = 0; r < mat_row; ++r) {
        row_elements += nzPerRow;
        srand(time(0));
        int *checkRepeat = (int *) calloc_or_exit(mat_row, sizeof(int));
        for (int i = 0; i < nzPerRow; i++) {
            int rand_idx;
            /// escape same random column
            do{
                rand_idx = rand() % mat_row;
            }while(checkRepeat[rand_idx]>0);
            checkRepeat[rand_idx] = 1;
            col_ptr[start_idx] = rand_idx;
            /// Fill by any random double value
            val_ptr[start_idx] = (double)(1 + (rand_idx %10));
            start_idx++;
        }
        row_ptr[r+1] = row_elements;
    }

    return 1;
}

/*
 * Tries to malloc. Terminates on failure.
 */
void * malloc_or_exit(size_t size) {
    void * ptr = malloc( size );
    if ( !ptr ) {
        fprintf(stderr, "malloc: failed to allocate memory\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/*
 * Tries to calloc. Terminates on failure.
 */
void * calloc_or_exit(size_t nmemb, size_t size) {
    void * ptr = calloc(nmemb, size);
    if ( !ptr ) {
        fprintf(stderr, "calloc: failed to allocate memory\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/*
 * Check if integer num is in range [start, start + count)
 */
char in_range(int num, int start, int count) {
    return (start <= num && num < start + count);
}

/*
 * Check if integer num is in range [start, start + count)
 */
char in_diagonal(int num, int start, int end) {
    return (start <= num && num <= end);
}


#include <time.h>

#include "util.h"

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
            int rand_idx;
            do{
                rand_idx = rand() % mat_row;
            }while(checkRepeat[rand_idx]>0);
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


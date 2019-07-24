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

int random_mat (int *buf_i, int *buf_j, double *buf_val, int start_row, int end_row, int nz)
{
    int num_row = end_row - start_row + 1;
    int nz_per_row = nz/num_row;
    int start_idx = 0;
    for (int r = 0; r < num_row; ++r) {
        int limit = nz_per_row-1;
        srand(time(0));
        for (int i = 0; i < nz_per_row; i++) {
            int idx = start_row + (rand() % num_row);
            buf_i[start_idx] = r;
            buf_j[start_idx] = idx;
            /// Fill by any random double value
            buf_val[start_idx] = (double)(1 + (idx %10));
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


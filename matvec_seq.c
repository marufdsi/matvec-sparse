#include <stdio.h>
#include <stdlib.h>

#define MAX_RANDOM_NUM (1<<20)

#include "mmio-wrapper.h"
#include "stopwatch.h"
#include "util.h"
#include <math.h>

#define TOTAL_RUNS 100

void mat_vec_mult(const double *values,
                    const int *i_idx,
                    const int *j_idx,
                    const double *x,
                    double *y,
                    const int NZ)
{
    for (int k = 0 ; k < NZ; k++) {
        y[ i_idx[k] ] += values[k] * x[ j_idx[k] ];
    }
}


int main(int argc, char * argv[])
{
    char *in_file,
         *out_file = NULL;

    double *values; /* a_values array */
    int *i_idx,     /* i_index array */
        *j_idx;     /* j_index array */

    double *x, *y;  /* Ax = y */
    double comp_time = 0;
    int N, NZ;

    /* read arguments */
    if (argc < 2 || argc > 3) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    }
    else {
        in_file = argv[1];
        if (argc == 3) 
            out_file = argv[2];
    }

    /* read matrix */
    if ( read_matrix(in_file, &i_idx, &j_idx, &values, &N, &NZ) != 0) {
        exit(EXIT_FAILURE);
    }

    debug("Matrix properties: N = %d, NZ = %d\n", N, NZ);

    /* allocate x, y vector */
    y = (double *)calloc( N, sizeof(double) );
    x = (double *)malloc( N * sizeof(double) );
    if (y == NULL || x == NULL) {
        fprintf(stderr, "malloc: failed to allocate memory\n");
        exit(EXIT_FAILURE);
    }

    /* generate random vector */
    //random_vec(x, N, MAX_RANDOM_NUM);
    for (int i = 0; i < N; i++) {
        x[i] = 1;
    }
    
    /* perform the multiplication */
    __sw_start(0);
    mat_vec_mult(values, i_idx, j_idx, x, y, NZ);
    __sw_stop(0, &comp_time);

    double stdev = 0, mean = 0, runs[TOTAL_RUNS];
    for (int r = 0; r < TOTAL_RUNS; r++) {
        double my_comp_time = 0;
        __sw_start(0);
        mat_vec_mult(values, i_idx, j_idx, x, y, NZ);
        //MPI_Barrier(MPI_COMM_WORLD);
        __sw_stop(0, &my_comp_time);
        runs[r] = (my_comp_time) * 1000.0;
        mean += runs[r];
    }
    printf("Total Runtime of %d times: %10.3lf\n", TOTAL_RUNS,  mean);
    mean /= TOTAL_RUNS;
    for (int r = 0; r < TOTAL_RUNS; r++) {
        stdev += (runs[r] - mean) * (runs[r] - mean);
    }
    stdev = sqrt(stdev);

    printf("Sequential Avg Computation time and Stdev: %10.3lf [%4.3lf] ms\n\n", mean, stdev);

    printf("Total execution time: %10.3lf ms\n", comp_time);

    if (out_file != NULL) {
        printf("Writing result to '%s'\n", out_file);

        /* open file */
        FILE *f;
        if ( !(f = fopen(out_file, "w")) ) {
            fprintf(stderr, "fopen: failed to open file '%s'", out_file);
            exit(EXIT_FAILURE);
        }

        /* write result */
        for (int i = 0; i < N; i++) {
            fprintf(f, "%.8lf\n", y[i]);
        }
        
        /* close file */
        if ( fclose(f) != 0) {
            fprintf(stderr, "fopen: failed to open file '%s'", out_file);
            exit(EXIT_FAILURE);
        }

        printf("Done!\n");
    }

    /* free the memory */
    free(values);
    free(i_idx);
    free(j_idx);
    free(x);
    free(y);
    
    return 0;
}

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "mmio-wrapper.h"
#include "util.h"
#include "partition.h"
#include <time.h>

void  create_random_matrix(int m, int n, int nnz_per_row, char *out_file){
    FILE *newMat;
    if (!(newMat = fopen(strcat(out_file, ".mtx"), "w"))) {
        printf("fopen: failed to open file '%s'", out_file);
        return;
    }
    fprintf(newMat, "%%%MatrixMarket matrix coordinate real general\n");
    fprintf(newMat, "%d %d %d\n", m, n, m*nnz_per_row);
    srand(time(NULL) * (n%m + 1));
    printf("*********** Initialization DOne **********");
    int *trackIndex;
    for (int i = 0; i < m; ++i) {
        trackIndex = (int *)calloc_or_exit( n, sizeof(int));
        for (int j = 0; j < nnz_per_row; ++j) {
            int randColIdx;
            do {
                randColIdx = rand() % n;
            } while (trackIndex[randColIdx] != 0);
            trackIndex[ randColIdx] = 1.0;
            fprintf(newMat, "%d %d %lf\n", i+1, randColIdx+1, ((double)(randColIdx%10) +1));
        }
        free(trackIndex);
    }
/// close file
    if (fclose(newMat) != 0) {
        fprintf(stderr, "fopen: failed to open file '%s'", out_file);
        return;
    }
}

void  create_diagonal_matrix(){

}

/**
 *
 * @param argc
 * @param argv
 * @return
 * Main method of the SpMV model
 * This model will able to predict the run time of the sparse matrix
 */
int main(int argc, char *argv[]) {

    char *out_file;
    int m, n, nnz_per_row;
    if (argc < 4) {
        printf("Usage: %s m n nnz_per_row output_file]\n", argv[0]);
        return 0;
    } else {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        nnz_per_row = atoi(argv[3]);
        out_file = argv[4];
    }

    create_random_matrix(m, n, nnz_per_row, out_file);

    return 0;
}


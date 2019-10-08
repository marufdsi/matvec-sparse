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

void  create_diagonal_matrix(int m, int n, int nnz_per_row, char *out_file){
    FILE *newMat;
    if (!(newMat = fopen(strcat(out_file, ".mtx"), "w"))) {
        printf("fopen: failed to open file '%s'", out_file);
        return;
    }
    fprintf(newMat, "%%%MatrixMarket matrix coordinate real general\n");
    fprintf(newMat, "%d %d %d\n", m, n, m*nnz_per_row);

    int start_idx = 0;
    int row_elements = 0;
    int lower_nnz = nnz_per_row - (nnz_per_row / 2);
    for (int r = 0; r < m; ++r) {
        int start_coldIdx = (r - lower_nnz + 1) < 0 ? 0 : ((r - lower_nnz + 1 + nnz_per_row)> m) ? (m-nnz_per_row) : (r - lower_nnz + 1);
        for (int colIdx = start_coldIdx; colIdx < start_coldIdx + nnz_per_row; ++colIdx) {
            fprintf(newMat, "%d %d %lf\n", r+1, colIdx+1, ((double)(colIdx%10) +1));
        }
    }
    /// close file
    if (fclose(newMat) != 0) {
        fprintf(stderr, "fopen: failed to open file '%s'", out_file);
        return;
    }
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
    int m, n, nnz_per_row, isDiagonal = 0;
    if (argc < 4) {
        printf("Usage: %s m n nnz_per_row output_file]\n", argv[0]);
        return 0;
    } else {
        m = atoi(argv[1]);
        n = atoi(argv[2]);
        nnz_per_row = atoi(argv[3]);
        out_file = argv[4];
        if(argc>5)
            isDiagonal = atoi(argv[5]);
    }
    if(isDiagonal == 0) {
        printf("Random Matrix Creation for %d rows, %d columns with %d non zeros per row.\n", m, n nnz_per_row);
        create_random_matrix(m, n, nnz_per_row, out_file);
    }
    else {
        printf("Diagonal Matrix Creation for %d rows, %d columns with %d non zeros per row.\n", m, n nnz_per_row);
        create_diagonal_matrix(m, n, nnz_per_row, out_file);
    }

    return 0;
}


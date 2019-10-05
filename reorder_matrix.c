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
#include <unordered_set>

int transactionByte = 128;

void
reOrderMatrix(int *row_ptr, int *col_ptr, ValueType *val_ptr, int m, int n, int nnzA, int max_deg, char *out_file) {
    FILE *newMat;
    if (!(newMat = fopen(strcat(out_file, ".mtx"), "w"))) {
        printf("fopen: failed to open file '%s'", out_file);
        return;
    }
    ValueType *x;
    x = (ValueType *) malloc(n * sizeof(ValueType));
    for (int i = 0; i < n; ++i) {
        x = 1;
    }
    fprintf(newMat, "%%%MatrixMarket matrix coordinate real general\n");
    fprintf(newMat, "%d %d %d\n", m, n, nnzA);
    int wordSize = transactionByte / sizeof(ValueType);
    for (int i = 0; i < m; ++i) {
        for (int deg = 0; deg < max_deg; ++deg) {
            int rows = (i + wordSize) > m ? m : (i + wordSize);
            std::unordered_set<long> hashme;
            int firstRow = i;
            while (firstRow < rows) {
                if (deg < (row_ptr[firstRow] - row_ptr[firstRow + 1])) {
                    break;
                }
                firstRow++;
            }
            if (firstRow == rows)
                break;
            hashme.insert((long) (&x[col_ptr[row_ptr[firstRow] + deg]]) / transactionByte);
            for (int w = firstRow + 1; w < rows; ++w) {
                if (deg < (row_ptr[w + 1] - row_ptr[w])) {
                    int match_count = 0, match_idx = row_ptr[w] + deg;
                    for (int pos = row_ptr[w] + deg; pos < row_ptr[w + 1]; ++pos) {
                        long x_add = (long) (&x[col_ptr[pos]]) / transactionByte;
                        if (hashme.count(x_add) > match_count) {
                            match_idx = col_ptr[pos];
                            match_count = hashme.count(x_add);
                        }
                    }
                    if (match_idx != row_ptr[w] + deg) {
                        int tmp = col_ptr[row_ptr[w] + deg];
                        col_ptr[row_ptr[w] + deg] = col_ptr[match_idx];
                        col_ptr[match_idx] = tmp;
                    }
                    hashme.insert((long) (&x[col_ptr[row_ptr[w] + deg]]) / transactionByte);
                }
            }
        }
    }

    for (int i = 0; i < m; ++i) {
        for (int adj = row_ptr[i]; adj < row_ptr[i + 1]; ++adj) {
            fprintf(newMat, "%d %d %lf\n", i, col_ptr[adj], val_ptr[adj]);
        }
    }
/// close file
    if (fclose(newMat) != 0) {
        fprintf(stderr, "fopen: failed to open file '%s'", out_file);
        return;
    }
}

void create_diagonal_matrix() {

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

    char *out_file, *in_file;
    int m, n, nnzA = 0, *row_ptr, *col_ptr;
    ValueType *val_ptr;
    if (argc < 2) {
        printf("Usage: %s input_file output_file]\n", argv[0]);
        return 0;
    } else {
        in_file = argv[1];
        if (argc > 2)
            out_file = argv[2];
    }

    int max_deg = 0;
    if (csr_read_Matrix(in_file, &row_ptr, &col_ptr, &val_ptr, &m, &n, &nnzA, &max_deg) != 0) {
        fprintf(stderr, "read_matrix: failed\n");
        exit(EXIT_FAILURE);
    }
    n = m;

    reOrderMatrix(row_ptr, col_ptr, val_ptr, m, n, nnzA, max_deg, out_file);

    return 0;
}


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "mmio-wrapper.h"

/**
 *
 * @param argc
 * @param argv
 * @return
 * Main method of the SpMV model
 * This model will able to predict the run time of the sparse matrix
 */
int main(int argc, char *argv[]) {

    char *in_file;
    int *row_ptr, *col_ptr, mat_row = 0, _nnz = 0;
    double *val_ptr;
    if (argc < 2) {
        printf("Usage: %s input_file [output_file]\n", argv[0]);
        return 0;
    } else {
        in_file = argv[1];
    }
    if (read_coo_matrix_to_csr(in_file, &row_ptr, &col_ptr, &val_ptr, &mat_row, &_nnz) != 0) {
        printf("Read Error for the file=%s\n", in_file);
    }

    if (mat_row <= 0) {
        printf("No rows are found for %s\n", in_file);
    }
    if (_nnz <= 0) {
        printf("No elements are found for %s\n", in_file);
    }

    printf("Find Mode\n");
    int maxValue = 0, maxCount = 0, max_deg=0, min_deg=INT_MAX;
    for (int i = 0; i < mat_row; ++i) {
        int deg = row_ptr[i+1]-row_ptr[i];
        if(deg>max_deg)
            max_deg = deg;
        if(deg<min_deg)
            min_deg = deg;

        int count = 0;
        for (int j = 0; j < mat_row; ++j) {
            if ((row_ptr[j + 1] - row_ptr[j]) == (row_ptr[i + 1] - row_ptr[i]))
                ++count;
        }
        if (count > maxCount) {
            maxCount = count;
            maxValue = row_ptr[i + 1] - row_ptr[i];
        }
    }

    char *ptr = strtok(in_file, "/");
    char *matrixName = strtok(strtok(NULL, "-"), ".");

    char outputFile[100] = "Matrix_Info.csv";
    FILE *resultCSV;
    FILE *checkFile;
    if ((checkFile = fopen(outputFile, "r")) != NULL) {
        // file exists
        fclose(checkFile);
        if (!(resultCSV = fopen(outputFile, "a"))) {
            fprintf(stderr, "fopen: failed to open %s file\n", outputFile);
            exit(EXIT_FAILURE);
        }
    } else {
        if (!(resultCSV = fopen(outputFile, "w"))) {
            fprintf(stderr, "fopen: failed to open file %s\n", outputFile);
            exit(EXIT_FAILURE);
        }
        fprintf(resultCSV, "MatrixName,MatrixSize,Mode_NNZ,Count,MinDegree,MaxDegree\n");
    }

    fprintf(resultCSV, "%s,%d,%d,%d,%d,%d\n", matrixName, mat_row, maxValue, maxCount, min_deg, max_deg);
    if (fclose(resultCSV) != 0) {
        fprintf(stderr, "fopen: failed to open file %s\n", outputFile);
        exit(EXIT_FAILURE);
    }

    return 0;
}


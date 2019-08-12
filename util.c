#include <time.h>
#include <jmorecfg.h>

#include "util.h"

double getVal(Map *map, int col, int limit) {
    if (limit > 0) {
        while (map != NULL && limit > 0) {
            if ((*map).key.col == col) {
                return (*map).value.val;
            }
            map++;
            limit--;
        }
    } else {
        while (map != NULL) {
            if ((*map).key.col == col) {
                return (*map).value.val;
            }
            map++;
        }
    }
//    printf("Error!!! column=%d not found\n", col);
    return -1.0;
}

/*
 * Generate a random real vector of size N with
 * elements beloning to [0, limit)
 */
void random_vec(double *v, int N, int limit) {
    /* fill v with random doubles */
    limit--;
    srand(410);
    for (int i = 0; i < N; i++) {
        v[i] = ((double) rand()) / (((double) RAND_MAX) / limit) + \
               ((double) rand()) / ((double) RAND_MAX);
    }
}

int random_mat(int *buf_i, int *buf_j, double *buf_val, int start_row, int mat_row, int nzPerRow) {
    int start_idx = 0;
    for (int r = 0; r < mat_row; ++r) {
        srand(time(0));
        int *checkRepeat = (int *) calloc_or_exit(mat_row, sizeof(int));
        for (int i = 0; i < nzPerRow; i++) {
//            int conflict = 0;
            int rand_idx;
            do {
//                conflict++;
                rand_idx = rand() % mat_row;
            } while (checkRepeat[rand_idx] > 0);
            /*if(conflict>0){
                printf("Conflicts Happened = %d\n", conflict);
            }*/
            checkRepeat[rand_idx] = 1;
            int idx = start_row + rand_idx;
            buf_i[start_idx] = start_row + r;
            buf_j[start_idx] = idx;
            /// Fill by any random double value
            buf_val[start_idx] = (double) (1 + (idx % 10));
            start_idx++;
        }
    }

    return 1;
}

int csr_random_mat(int rank, int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int mat_col, int nzPerRow,
                   int dist, int *non_diagonal_elements) {
    int col_idx = 0;
    int row_idx = 0;
    row_ptr[0] = row_idx;
    (*non_diagonal_elements) = 0;
    for (int r = 0; r < mat_row; ++r) {
        int *trackIndex = (int *) malloc_or_exit(mat_row * sizeof(int));
        int *isIdTaken = (int *) calloc_or_exit(mat_col, sizeof(int));
        for (int k = 0; k < mat_row; ++k) {
            trackIndex[k] = k;
        }
        int idx_not_taken = mat_row;
        row_idx += nzPerRow;
        int non_diagonal = 0;
        int id_range = mat_col;
        int start_id = 0;
        srand(time(NULL) * (r + 1) * (rank + 1));
        for (int i = 0; i < nzPerRow; i++) {
            if ((nzPerRow * dist) / 100 <= non_diagonal) {
                id_range = mat_row;
                start_id = rank * mat_row;
            }
            int rand_id;
            int maxTry = 0;
            /// escape same random column
            do {
                rand_id = rand() % id_range;
                maxTry++;
            } while (isIdTaken[rand_id] != 0 && maxTry < 5);
            if (maxTry >= 5 && isIdTaken[rand_id] != 0) {
                if (idx_not_taken <= 0) {
                    printf("[%d] Exception occurred\n", rank);
                    return 0;
                }
                int _idx = rand() % idx_not_taken;
                rand_id = trackIndex[_idx];
                trackIndex[_idx] = trackIndex[idx_not_taken - 1];
                idx_not_taken--;
            } else {
                if (rand_id < mat_row) {
                    if (rand_id == trackIndex[rand_id]) {
                        trackIndex[rand_id] = trackIndex[idx_not_taken - 1];
                        idx_not_taken--;
                    } else {
                        for (int l = 0; l < idx_not_taken; ++l) {
                            if (rand_id == trackIndex[l]) {
                                trackIndex[l] = trackIndex[idx_not_taken - 1];
                                idx_not_taken--;
                                break;
                            }
                        }
                    }
                }
            }
            isIdTaken[rand_id] = 1;
            if (!in_diagonal((start_id + rand_id), rank * mat_row, ((rank+1) * mat_row) - 1)) {
                non_diagonal++;
            }
            col_ptr[col_idx] = start_id + rand_id;
            /// Fill by any random double value
            val_ptr[col_idx] = (double) (1 + (rand_id % 10));
            col_idx++;
        }
        (*non_diagonal_elements) += non_diagonal;
        free(trackIndex);
        free(isIdTaken);
        row_ptr[r + 1] = row_idx;
    }

    return 1;
}

int csr_random_diagonal_mat(int rank, int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int nzPerRow) {
    int start_idx = 0;
    int row_elements = 0;
    row_ptr[0] = row_elements;
    int range_start = rank * mat_row;
    for (int r = 0; r < mat_row; ++r) {
        int *trackIndex = (int *) calloc_or_exit(mat_row, sizeof(int));
        row_elements += nzPerRow;
        Map *map = (Map *) malloc_or_exit(nzPerRow * sizeof(Map));
        srand(time(NULL) * (r + 1) * (rank + 1));
        int taken_idx = 0;
        for (int i = 0; i < nzPerRow; i++) {
            int rand_idx;
            int maxTry = 0;
            /// escape same random column
            do {
                maxTry++;
                rand_idx = rand() % mat_row;
            } while (!(getVal(map, rand_idx, i + 1) < 0) && maxTry < 50);
            if (maxTry >= 50) {
                for (int k = taken_idx + 1; k < mat_row; ++k) {
                    if (trackIndex[k] == 0) {
                        rand_idx = k;
                        taken_idx = k;
                        break;
                    }
                }
            }
            trackIndex[rand_idx] = 1;
            map[i].key.col = rand_idx;
            map[i].value.val = 1.0;
            col_ptr[start_idx] = range_start + rand_idx;
            /// Fill by any random double value
            val_ptr[start_idx] = (double) (1 + (rand_idx % 10));
            start_idx++;
        }
        free(trackIndex);
        free(map);
        row_ptr[r + 1] = row_elements;
    }

    return 1;
}

int csr_diagonal_mat(int rank, int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int nzPerRow) {
    int start_idx = 0;
    int row_elements = 0;
    row_ptr[0] = row_elements;
    int lower_nnz = nzPerRow - (nzPerRow / 2);
    int range_start = rank * mat_row;
    for (int r = 0; r < mat_row; ++r) {
        row_elements += nzPerRow;
        int start_coldIdx = (r - lower_nnz + 1) < 0 ? 0 : ((r - lower_nnz + 1 + nzPerRow)> mat_row) ? (mat_row-nzPerRow) : (r - lower_nnz + 1);
        for (int colIdx = start_coldIdx; colIdx < nzPerRow + start_coldIdx; ++colIdx) {
            col_ptr[start_idx] = range_start + colIdx;
            /// Fill by any random double value
            val_ptr[start_idx] = (double) (1 + (colIdx % 10));
            start_idx++;
        }
        row_ptr[r + 1] = row_elements;
    }
    return 1;
}

int csr_diagonal_mat_with_bandwidth(int rank, int *row_ptr, int *col_ptr, double *val_ptr, int mat_row, int nzPerRow, int gap, int *bandwidth) {
    int start_idx = 0;
    int row_elements = 0;
    row_ptr[0] = row_elements;
    (*bandwidth) = nzPerRow + (nzPerRow-1)*gap;
    int lower_nnz = (*bandwidth) - ((*bandwidth) / 2);
    int range_start = rank * mat_row;
    for (int r = 0; r < mat_row; ++r) {
        row_elements += nzPerRow;
        int start_coldIdx = (r - lower_nnz + 1) < 0 ? 0 : ((r - lower_nnz + 1 + (*bandwidth))> mat_row) ? (mat_row-(*bandwidth)) : (r - lower_nnz + 1);
        for (int colIdx = start_coldIdx; colIdx < (*bandwidth) + start_coldIdx; colIdx+=gap) {
            col_ptr[start_idx] = range_start + colIdx;
            /// Fill by any random double value
            val_ptr[start_idx] = (double) (1 + (colIdx % 10));
            start_idx ++;
        }
        row_ptr[r + 1] = row_elements;
    }
    return 1;
}


/*
 * Tries to malloc. Terminates on failure.
 */
void *malloc_or_exit(size_t size) {
    void *ptr = malloc(size);
    if (!ptr) {
        fprintf(stderr, "malloc: failed to allocate memory\n");
        exit(EXIT_FAILURE);
    }
    return ptr;
}

/*
 * Tries to calloc. Terminates on failure.
 */
void *calloc_or_exit(size_t nmemb, size_t size) {
    void *ptr = calloc(nmemb, size);
    if (!ptr) {
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


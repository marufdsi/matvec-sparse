#!/bin/bash

export OMP_NUM_THREADS=$thread
export KMP_AFFINITY=granularity=fine,balanced,0,0

./omp_spmv_model $row $nnz $run $type "balanced" 0

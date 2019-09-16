#!/bin/bash


export OMP_NUM_THREADS=$thread
export KMP_AFFINITY=granularity=fine,scatter,0,0


./omp_spmv_model $row $nnz $run $type "scatter" 0

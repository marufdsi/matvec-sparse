#!/bin/bash

export OMP_NUM_THREADS=$thread
export KMP_AFFINITY=granularity=fine,$affinity,0,0


./omp_spmv input/nlpkkt160_original.mtx $run $affinity $reserve
./omp_spmv input/nlpkkt240_original.mtx $run $affinity $reserve
./omp_spmv input/cage15_original.mtx $run $affinity $reserve
./omp_spmv input/G3_circuit_original.mtx $run $affinity $reserve
./omp_spmv input/kkt_power_original.mtx $run $affinity $reserve
./omp_spmv input/ecology1_original.mtx $run $affinity $reserve
./omp_spmv input/uk-2002_original.mtx $run $affinity $reserve
./omp_spmv input/asia_original.mtx $run $affinity $reserve
./omp_spmv input/NLR_original.mtx $run $affinity $reserve
./omp_spmv input/coPapersDBLP_original.mtx $run $affinity $reserve
./omp_spmv input/coPapersCiteseer_original.mtx $run $affinity $reserve
./omp_spmv input/copter2_original.mtx $run $affinity $reserve
./omp_spmv input/333SP_original.mtx $run $affinity $reserve
./omp_spmv input/af_shell10_original.mtx $run $affinity $reserve
./omp_spmv input/AS365_original.mtx $run $affinity $reserve
./omp_spmv input/europe_original.mtx $run $affinity $reserve
./omp_spmv input/M6_original.mtx $run $affinity $reserve
./omp_spmv input/nlpkkt200_original.mtx $run $affinity $reserve
./omp_spmv input/thermal2_original.mtx $run $affinity $reserve

#./omp_spmv $file $run $affinity $reserve

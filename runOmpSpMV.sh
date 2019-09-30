#!/bin/bash

export OMP_NUM_THREADS=$thread
export KMP_AFFINITY=granularity=fine,$affinity,0,0


#./omp_spmv input/nlpkkt160_900.mtx $run $affinity $reserve
#./omp_spmv input/nlpkkt240_900.mtx $run $affinity $reserve
#./omp_spmv input/cage15_900.mtx $run $affinity $reserve
#./omp_spmv input/G3_circuit_900.mtx $run $affinity $reserve
#./omp_spmv input/kkt_power_900.mtx $run $affinity $reserve
#./omp_spmv input/ecology1_900.mtx $run $affinity $reserve
#./omp_spmv input/uk-2002_900.mtx $run $affinity $reserve
#./omp_spmv input/asia_900.mtx $run $affinity $reserve
#./omp_spmv input/NLR_900.mtx $run $affinity $reserve
#./omp_spmv input/coPapersDBLP_900.mtx $run $affinity $reserve
#./omp_spmv input/coPapersCiteseer_900.mtx $run $affinity $reserve
#./omp_spmv input/copter2_900.mtx $run $affinity $reserve
#./omp_spmv input/333SP_900.mtx $run $affinity $reserve
#./omp_spmv input/af_shell10_900.mtx $run $affinity $reserve
#./omp_spmv input/AS365_900.mtx $run $affinity $reserve
#./omp_spmv input/europe_900.mtx $run $affinity $reserve
#./omp_spmv input/M6_900.mtx $run $affinity $reserve
#./omp_spmv input/nlpkkt200_900.mtx $run $affinity $reserve
#./omp_spmv input/thermal2_900.mtx $run $affinity $reserve

./omp_spmv $file $run $affinity $reserve

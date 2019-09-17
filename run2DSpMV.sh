#!/bin/bash
module load openmpi

mpirun -n $procs ./spmv_random "input/"$input"_random_"$procs".mtx"


#mpirun -n $procs ./spmv_random "input/nlpkkt160_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/nlpkkt240_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/cage15_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/G3_circuit_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/kkt_power_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/ecology1_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/uk-2002_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/asia_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/NLR_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/coPapersDBLP_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/coPapersCiteseer_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/copter2_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/333SP_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/af_shell10_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/AS365_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/europe_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/M6_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/nlpkkt200_random_"$procs".mtx"
#mpirun -n $procs ./spmv_random "input/thermal2_random_"$procs".mtx"

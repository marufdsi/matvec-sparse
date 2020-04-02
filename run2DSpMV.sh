#!/bin/bash
module remove openmpi/4.0.1
module load openmpi/3.1.2

#mpirun -n $procs ./spmv_random "input/"$input"_random_"$procs".mtx" $nodes

#mpirun -n $procs ./spmv_random "input/cnr_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/eu_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/in_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/NACA0015_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/road_usa_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/road_central_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/nlpkkt120_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/nlpkkt160_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/nlpkkt200_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/nlpkkt240_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/cage15_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/G3_circuit_random_"$procs".mtx" $nodes
mpirun -n $procs ./spmv_random "input/kkt_power_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/ecology1_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/uk2002_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/uk2007_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/asia_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/NLR_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/coPapersDBLP_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/citationCiteseer_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/coPapersCiteseer_random_"$procs".mtx" $nodes
##mpirun -n $procs ./spmv_random "input/copter2_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/333SP_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/af_shell10_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/AS365_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/europe_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/M6_random_"$procs".mtx" $nodes
#mpirun -n $procs ./spmv_random "input/thermal2_random_"$procs".mtx" $nodes

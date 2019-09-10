#!/bin/bash
module load openmpi

#mpirun -n $procs ./spmv_p2p "input/"$input"_"$procs".mtx"


mpirun -n $procs ./mult_p2p "input/nlpkkt160_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/nlpkkt240_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/cage15_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/G3_circuit_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/kkt_power_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/ecology1_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/uk-2002_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/asia_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/NLR_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/coPapersDBLP_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/coPapersCiteseer_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/copter2_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/333SP_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/af_shell10_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/AS365_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/europe_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/M6_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/nlpkkt200_"$procs".mtx"
mpirun -n $procs ./mult_p2p "input/thermal2_"$procs".mtx"


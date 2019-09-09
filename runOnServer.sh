#!/bin/bash
module load openmpi
#module load mpich
echo "input/"$graph"_"$parts".mtx"
#mpirun -np $thread ./matvec_mpi_p2p "input/"$graph"_"$parts"_original.mtx" "output/"$graph"_"$parts"_original.mtx"

mpirun -n $thread ./matvec_mpi_p2p "input/"$graph"_"$parts".mtx" "output/"$graph"_"$parts".mtx"

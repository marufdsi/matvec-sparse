#!/bin/bash
module load openmpi

mpirun -n $thread ./matvec_mpi_calculation $size $nz $run

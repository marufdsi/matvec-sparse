#!/bin/bash
module load openmpi

mpirun -n $thread ./csr_mpi_model $row $col $nz $run

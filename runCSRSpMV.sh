#!/bin/bash
module load openmpi

mpirun -n $thread ./csr_mpi_spmv $row $col $nz $run



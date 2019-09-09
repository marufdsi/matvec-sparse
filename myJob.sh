#!/bin/bash

echo "************** Original Matrix *********************"

mpirun -np 2 ./matvec_mpi_p2p input/ecology1_2_original.mtx output/ecology1_2_original.mtx



echo "******************* Partitioned Matrix ********************"

#mpirun -np 2 ./matvec_mpi_p2p input/ecology1_2.mtx output/ecology1_2.mtx

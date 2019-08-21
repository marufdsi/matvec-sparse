cmake_minimum_required(VERSION 3.9)
project(matvec_sparse C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(matvec_sparse
        matvec_mpi_bcast.c
        matvec_mpi_calculation.c
        matvec_mpi_p2p.c
        matvec_seq.c
        mmio-wrapper.c
        mmio-wrapper.h
        mmio.c
        mmio.h
        partition.c
        partition.h
        stopwatch.c
        stopwatch.h
        util.c
        util.h
        csr_mpi_spmv.c
        csr_mpi_model.c
        csr_random_spmv_model.c
        csr_spmv.c
        csr_mpi_reduced_spmv.c
        spmv_p2p.c
        comm_p2p.c)

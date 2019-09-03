# C compiler 
CC = gcc
CFLAGS = -O3 -std=gnu99 -DDEBUG

OMP_FLAGS = -fopenmp

# MPI compiler wrapper
MPI_C = mpicc
MPI_CFLAGS = -O3 -std=gnu99 -DDEBUG

# Link libraries
LDFLAGS = -lm -lrt

# Object files from libraries
OBJ = mmio.o mmio-wrapper.o partition.o util.o

all: matvec_seq matvec_mpi_p2p matvec_mpi_bcast matvec_mpi_calculation csr_mpi_spmv csr_mpi_model csr_random_spmv_model csr_spmv csr_mpi_reduced_spmv spmv_p2p comm_p2p mult_p2p spmv_random read_file spmv_random_model omp_spmv_model

matvec_seq: matvec_seq.c $(OBJ) stopwatch.o
	$(CC) $(CFLAGS) $(OBJ) stopwatch.o $< -o $@ $(LDFLAGS)

matvec_mpi_p2p: matvec_mpi_p2p.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

matvec_mpi_bcast: matvec_mpi_bcast.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

matvec_mpi_calculation: matvec_mpi_calculation.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

csr_mpi_spmv: csr_mpi_spmv.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

csr_mpi_model: csr_mpi_model.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

csr_random_spmv_model: csr_random_spmv_model.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

csr_spmv: csr_spmv.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

csr_mpi_reduced_spmv: csr_mpi_reduced_spmv.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

spmv_p2p: spmv_p2p.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

comm_p2p: comm_p2p.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

mult_p2p: mult_p2p.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

spmv_random: spmv_random.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

spmv_random_model: spmv_random_model.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

read_file: read_file.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

omp_spmv_model: omp_spmv_model.c $(OBJ)
	$(CC) $(CFLAGS) $(OMP_FLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f matvec_seq matvec_mpi_bcast matvec_mpi_p2p matvec_mpi_calculation csr_mpi_spmv csr_mpi_model csr_random_spmv_model csr_spmv csr_mpi_reduced_spmv spmv_p2p comm_p2p mult_p2p spmv_random read_file spmv_random_model omp_spmv_model *.o

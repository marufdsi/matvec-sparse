# C compiler 
CC = gcc
CFLAGS = -O3 -std=gnu99 -DDEBUG

# MPI compiler wrapper
MPI_C = mpicc
MPI_CFLAGS = -O3 -std=gnu99 -DDEBUG

# Link libraries
LDFLAGS = -lm -lrt

# Object files from libraries
OBJ = mmio.o mmio-wrapper.o partition.o util.o

all: matvec_seq matvec_mpi_p2p matvec_mpi_bcast matvec_mpi_calculation csr_mpi_spmv csr_mpi_model csr_spmv

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

csr_spmv: csr_spmv.c $(OBJ)
	$(MPI_C) $(MPI_CFLAGS) $(OBJ) $< -o $@ $(LDFLAGS)

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f matvec_seq matvec_mpi_bcast matvec_mpi_p2p matvec_mpi_calculation csr_mpi_spmv csr_mpi_model csr_spmv *.o

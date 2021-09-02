CC=gcc
MPICC=mpicc
#CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3 -std=c99 

#default: all
# OMP implementations
main_simple3: mmio.o main_simple3.o 
	$(CC) $(CFLAGS) -o $@ mmio.o main_simple3.o -fopenmp

main_simple3.o: main_simple3.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_simple3.c -fopenmp

# MPI Implementation
main_mpi: mmio.o main_mpi.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpi.o

main_mpi.o: main_mpi.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpi.c

# MPI + OMP Hybrid Implementation
main_hyb: mmio.o main_hyb.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_hyb.o -fopenmp

main_hyb.o: main_hyb.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_hyb.c -fopenmp
#---------------------------------------------
# Compiling libraries
mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

.PHONY: clean

clean:
	rm -f main_simple3 main_simple3.o mmio.o
	rm -f main_mpi main_mpi.o mmio.o
	rm -f main_hyb main_hyb.o mmio.o

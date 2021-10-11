CC=gcc
MPICC=mpicc
#CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3 -std=c99 

#default: all
# OMP implementations
main_omp: mmio.o main_omp.o 
	$(CC) $(CFLAGS) -o $@ mmio.o main_omp.o -fopenmp

main_omp.o: main_omp.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_omp.c -fopenmp

main_omp2: mmio.o main_omp2.o 
	$(CC) $(CFLAGS) -o $@ mmio.o main_omp2.o -fopenmp

main_omp2.o: main_omp2.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_omp2.c -fopenmp
# MPI Implementation
main_mpi: mmio.o main_mpi.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpi.o

main_mpi.o: main_mpi.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpi.c

# MPI Implementation Version2
main_mpi2: mmio.o main_mpi2.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpi2.o

main_mpi2.o: main_mpi2.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpi2.c


# MPI + OMP Hybrid Implementation
main_hyb: mmio.o main_hyb.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_hyb.o -fopenmp

main_hyb.o: main_hyb.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_hyb.c -fopenmp

main_hyb2: mmio.o main_hyb2.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_hyb2.o -fopenmp

main_hyb2.o: main_hyb2.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_hyb2.c -fopenmp
#---------------------------------------------
# Compiling libraries
mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

.PHONY: clean

clean:
	rm -f main_omp main_omp.o mmio.o
	rm -f main_omp2 main_omp2.o mmio.o
	rm -f main_mpi main_mpi.o mmio.o
	rm -f main_mpi2 main_mpi2.o mmio.o
	rm -f main_hyb main_hyb.o mmio.o
	rm -f main_hyb2 main_hyb2.o mmio.o

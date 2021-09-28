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

main_omp: mmio.o main_omp.o 
	$(CC) $(CFLAGS) -o $@ mmio.o main_omp.o -fopenmp

main_omp.o: main_omp.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_omp.c -fopenmp

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

# MPI Implementation Version3
main_mpi3: mmio.o main_mpi3.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpi3.o

main_mpi3.o: main_mpi3.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpi3.c


# MPI Implementation Version4
main_mpi4: mmio.o main_mpi4.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpi4.o

main_mpi4.o: main_mpi4.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpi4.c

# MPI Implementation Version4_1
main_mpi4_1: mmio.o main_mpi4_1.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpi4_1.o

main_mpi4_1.o: main_mpi4_1.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpi4_1.c

# MPI Implementation Version5
main_mpi5: mmio.o main_mpi5.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpi5.o

main_mpi5.o: main_mpi5.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpi5.c
# MPI + OMP Hybrid Implementation
main_hyb: mmio.o main_hyb.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_hyb.o -fopenmp

main_hyb.o: main_hyb.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_hyb.c -fopenmp

main_hyb_1: mmio.o main_hyb_1.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_hyb_1.o -fopenmp

main_hyb_1.o: main_hyb_1.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_hyb_1.c -fopenmp

#---------------------------------------------
# Compiling libraries
mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

.PHONY: clean

clean:
	rm -f main_simple3 main_simple3.o mmio.o
	rm -f main_omp main_omp.o mmio.o
	rm -f main_mpi main_mpi.o mmio.o
	rm -f main_mpi2 main_mpi2.o mmio.o
	rm -f main_mpi3 main_mpi3.o mmio.o
	rm -f main_mpi4 main_mpi4.o mmio.o
	rm -f main_mpi4_1 main_mpi4_1.o mmio.o
	rm -f main_mpi5 main_mpi5.o mmio.o
	rm -f main_hyb main_hyb.o mmio.o
	rm -f main_hyb_1 main_hyb_1.o mmio.o

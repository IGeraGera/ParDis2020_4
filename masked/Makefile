CC=gcc
MPICC=mpicc
#CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-O3 -std=c99 

#default: all
# OMP implementations

main_ompfinal: mmio.o main_ompfinal.o 
	$(CC) $(CFLAGS) -o $@ mmio.o main_ompfinal.o -fopenmp

main_ompfinal.o: main_ompfinal.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_ompfinal.c -fopenmp

# MPI Implementation Version
main_mpifinal: mmio.o main_mpifinal.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpifinal.o

main_mpifinal.o: main_mpifinal.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpifinal.c


main_mpiblocks: mmio.o main_mpiblocks.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_mpiblocks.o

main_mpiblocks.o: main_mpiblocks.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_mpiblocks.c

# MPI + OMP Hybrid Implementation
main_hybfinal: mmio.o main_hybfinal.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_hybfinal.o -fopenmp

main_hybfinal.o: main_hybfinal.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_hybfinal.c -fopenmp

main_hybblocks: mmio.o main_hybblocks.o
	$(MPICC) $(CFLAGS) -o $@ mmio.o main_hybblocks.o -fopenmp

main_hybblocks.o: main_hybblocks.c inc/mmio.h
	$(MPICC) $(CFLAGS) -c main_hybblocks.c -fopenmp

#---------------------------------------------
# Compiling libraries
mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

.PHONY: clean

clean:
	rm -f main_ompfinal main_ompfinal.o mmio.o
	rm -f main_mpifinal main_mpifinal.o mmio.o
	rm -f main_hybfinal main_hybfinal.o mmio.o
	rm -f main_mpiblocks main_mpiblocks.o mmio.o
	rm -f main_hybblocks main_hybblocks.o mmio.o

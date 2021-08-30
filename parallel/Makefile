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

#---------------------------------------------
# Compiling libraries
mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

.PHONY: clean

clean:
	rm -f main_simple3 main_simple3.o mmio.o

CC=gcc
CFLAGS=-O3 -std=c99

main: mmio.o main.o
	$(CC) $(CFLAGS) -o $@ mmio.o main.o

main.o: main.c inc/mmio.h
	$(CC) $(CFLAGS) -c main.c


mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

.PHONY: clean

clean:
	-rm -f mmio.o
	-rm -f main.o
	-rm -f main


CC=gcc
CFLAGS=-O3 -std=c99

all: main main_simple

main: mmio.o main.o
	$(CC) $(CFLAGS) -o $@ mmio.o main.o

main.o: main.c inc/mmio.h
	$(CC) $(CFLAGS) -c main.c

main_simple: mmio.o main_simple.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_simple.o

main_simple.o: main_simple.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_simple.c

main_simple2: mmio.o main_simple2.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_simple2.o

main_simple2.o: main_simple2.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_simple2.c

main_simpleCSR: mmio.o main_simpleCSR.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_simpleCSR.o

main_simpleCSR.o: main_simpleCSR.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_simpleCSR.c

mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

.PHONY: clean

clean:
	-rm -f mmio.o
	-rm -f main.o
	-rm -f main
	-rm -f main_simple.o
	-rm -f main_simple
	-rm -f main_simple2.o
	-rm -f main_simple2
	-rm -f main_simpleCSR.o
	-rm -f main_simpleCSR


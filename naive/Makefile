CC=gcc
CFLAGS=-O3 -std=c99
DBFLAGS=-O0 -std=c99 -g

all: main main_simple

debug: main_4R_LUT_db 

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

main_simple3: mmio.o main_simple3.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_simple3.o

main_simple3.o: main_simple3.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_simple3.c

main_simpleCSR: mmio.o main_simpleCSR.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_simpleCSR.o

main_simpleCSR.o: main_simpleCSR.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_simpleCSR.c

main_4R_LUT: mmio.o main_4R_LUT.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_4R_LUT.o 

main_4R_LUT.o: main_4R_LUT.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_4R_LUT.c 

main_naive: mmio.o main_naive.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_naive.o

main_naive.o: main_naive.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_naive.c

mmio.o: inc/mmio.c
	$(CC) $(CFLAGS) -c inc/mmio.c

# Debug Area

main_4R_LUT_db: mmio.o main_4R_LUT_db.o
	$(CC) $(DBFLAGS) -o $@ mmio.o main_4R_LUT_db.o 

main_4R_LUT_db.o: main_4R_LUT.c inc/mmio.h
	$(CC) $(DBLAGS) -o $@ -c main_4R_LUT.c 

.PHONY: clean

clean:
	-rm -f mmio.o
	-rm -f main.o
	-rm -f main
	-rm -f main_simple.o
	-rm -f main_simple
	-rm -f main_naive.o
	-rm -f main_naive
	-rm -f main_4R_LUT.o
	-rm -f main_4R_LUT
	-rm -f main_simple2.o
	-rm -f main_simple2
	-rm -f main_simple3.o
	-rm -f main_simple3
	-rm -f main_simpleCSR.o
	-rm -f main_simpleCSR
	-rm -f main_4R_LUT_db.o
	-rm -f main_4R_LUT_db



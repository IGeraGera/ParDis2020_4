CC=gcc
CFLAGS=-O3 -std=c99

all: main main_simple



main_naive: mmio.o main_naive.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_naive.o

main_naive.o: main_naive.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_naive.c

main_naive_masked: mmio.o main_naive_masked.o
	$(CC) $(CFLAGS) -o $@ mmio.o main_naive_masked.o

main_naive_masked.o: main_naive_masked.c inc/mmio.h
	$(CC) $(CFLAGS) -c main_naive_masked.c

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
	-rm -f 	-rm -f main_naive main_naive.o
	-rm -f main_naive_masked.o main_naive_masked



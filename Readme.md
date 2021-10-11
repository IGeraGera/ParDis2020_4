# Parallel & Distibuted Systems Assignment 4

## Info

Naive contains the 2 naive implementations for C=A\*B and C=F.\*(A\*B)


Parallel contains the MPI, OpenMP and Hybrid implementations for the final method and the blocks method for C=A\*B 

masked contains the MPI, OpenMP and Hybrid implementations for the final method and the blocks method for C=F.\*(A\*B)

## How to compile

Dependencies: 

+ mpicc (for MPI)
+ fopenmp (for OpenMP)

For example:

``` bash
cd naive
make main_naive
./main_naive ../dataset/test_ASmall.mtx ../dataset/test_BSmall.mtx
make clean main_naive_masked
./main_naive_masked ../dataset/test_ASmall.mtx ../dataset/test_BSmall.mtx ../dataset/test_FSmall.mtx
```
 


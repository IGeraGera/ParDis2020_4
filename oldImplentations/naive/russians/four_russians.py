import numpy as np
import math as mt
from numpy.core.shape_base import hstack

file1 = 'array_1.txt'
file2 = 'array_2.txt'
file3 = 'array_3.txt'

A = np.loadtxt(file1,dtype=bool,delimiter=',')
B = np.loadtxt(file2,dtype=bool,delimiter=',')
C_correct = np.loadtxt(file3,dtype=bool,delimiter=',')
C = np.zeros_like(A)

total = 13
b = 3
nb = mt.ceil(13/3)

for i in range(nb):
    row = i * b
    row_end =  row + b
    if row_end> total:
        Ai=A[:,row:total]
        Ai= np.hstack((Ai,np.zeros([total,row_end-total],dtype=bool)))    
    else:
        Ai=A[:,row:row_end]
    col = i * b
    col_end = col + b
    if col_end>total:
        Bi=B[col:total,:]
        Bi= np.vstack((Bi,np.zeros([col_end-total,total],dtype=bool)))
        print(Bi.shape)
    else:
        Bi=B[col:col_end,:]
    C = np.logical_or(np.matmul(Ai,Bi), C)

print(C == C_correct)


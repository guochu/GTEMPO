from triqs.gf import *
from triqs.operators import *
from h5 import *
from triqs_cthyb import Solver

A = HDFArchive("bethe.h5", 'r')

for i in range(0,10):
    f = open("Giw-%i.dat"%(i+1), 'w')
    G = A['G-%i'%i]
    for iw in range (0,2050):
        up = G['up'].data[iw,0,0]
        print(2*(iw-1024)-1, up.real, up.imag, file=f)
    f.close()

    
for i in range(0,10):
    f = open("Gtau-%i.dat"%(i+1), 'w')
    G = A['G_tau-%i'%i]
    for t in range(0,10001):
        up = G['up'].data[t,0,0]
        print(t/1000, up.real, up.imag, file=f)
    f.close()
        

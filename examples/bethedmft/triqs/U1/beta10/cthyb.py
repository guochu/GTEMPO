from triqs.gf import *
from triqs.operators import *
from h5 import *
import triqs.utility.mpi as mpi
from triqs_cthyb import Solver

# Parameters of the model
U = 1.0
t = 0.5
mu = U/2.0
beta = 10.0
n_loops = 10

p = {}
p["max_time"] = -1
p["length_cycle"] = 10
p["n_warmup_cycles"] = 10000
p["n_cycles"] = 10**8
p["use_norm_as_weight"] = True
p["measure_density_matrix"] = False
p["move_double"] = True

# Construct the impurity solver
S = Solver(beta = beta, gf_struct = [('up',1), ('down',1)] )

# This is a first guess for G
S.G_iw << t*SemiCircular(2*t)

# DMFT loop with self-consistency
for i in range(n_loops):

    if mpi.is_master_node():
        print("\n\nIteration = %i / %i" % (i+1, n_loops))

    # Symmetrize the Green's function imposing paramagnetism and use self-consistency
    g = 0.5 * ( S.G_iw['up'] + S.G_iw['down'] )
    for name, g0 in S.G0_iw:
        g0 << inverse( iOmega_n + mu - t**2 * g )

    # Solve the impurity problem
    S.solve(h_int = U * n('up',0) * n('down',0), **p)

    # Save iteration in archive
    if mpi.is_master_node():
       with HDFArchive("bethe.h5", 'a') as A:
           A['G_tau-%i'%i] = S.G_tau
           A['G-%i'%i] = S.G_iw
           A['Sigma-%i'%i] = S.Sigma_iw

    #if mpi.is_master_node():
        # f = open("Giw-%i.dat"%i, 'w')
        # A = S.G_iw
        # for iw in range (0,2050):
        #     up = A['up'].data[iw,0,0]
        #     print(2*(iw-1024)-1, up.real, up.imag, file=f)
        # f.close()
        # f = open("Gtau-%i.dat"%i, 'w')
        # A = S.G_tau
        # for t in range(0,10001):
        #     up = A['up'].data[t,0,0]
        #     print(t/10000.0, up.real, up.imag, file=f)
        # f.close()
        # f = open("Sigma-%i.dat"%i, 'w')
        # A = S.Sigma_iw
        # for iw in range (0,2050):
        #     up = A['up'].data[iw,0,0]
        #     print(2*(iw-1024)-1, up.real, up.imag, file=f)
        # f.close()

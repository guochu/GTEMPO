import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math


rc('text', usetex=True)


def read_tempo(beta,i):
	filename = 'beta%s/Giw-%s.dat'%(beta, i)
	data = loadtxt(filename)
	return data[:, 1]


def read_mc(beta, i):
	mpath = '../triqs/U1/beta%s/'%(beta)
	filename = mpath + 'Giw-%s.dat'%(i)
	data = loadtxt(filename)
	return data[:, 2]

fontsize = 20
labelsize = 16
linewidth = 2.
markersize = 10

colors = ['b', 'c', 'orange', 'g', 'y', 'r', 'k']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'r'

beta = 10

data1 = [read_tempo(beta, step) for step in range(1, 11)]
data2 = [read_mc(beta, step) for step in range(1, 11)]

distances = [norm(a-b)/norm(a) for (a, b) in zip(data1, data2)]

fig, ax = plt.subplots(1,1, figsize=(8, 7))


ax.plot(range(1, 11), distances, ls='--', marker='o', markersize=10, markerfacecolor='none', color='r', linewidth=2)



plt.tight_layout(pad=0.5)

# # plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


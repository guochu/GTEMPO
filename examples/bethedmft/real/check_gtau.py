import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math


rc('text', usetex=True)


def read_tempo(beta, i):
	filename = 'beta%s/Gtau-%s.dat'%(beta, i)
	data = loadtxt(filename)
	return data[:, 1]


def read_mc(beta, i):
	mpath = '../triqs/U1/beta%s/'%(beta)
	filename = mpath + 'Gtau-%s.dat'%(i)
	data = loadtxt(filename)
	return data[:, 1]

fontsize = 20
labelsize = 16
linewidth = 2.
markersize = 10

colors = ['b', 'c', 'orange', 'g', 'y', 'r', 'k']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'r'


beta = 1
step = 1

data1 = read_tempo(beta, step)
betas1 = linspace(0, beta, num=len(data1))
data2 = read_mc(beta, step)
betas2 = linspace(0, beta, num=len(data2))

fig, ax = plt.subplots(1,1, figsize=(8, 7))


ax.plot(betas1, data1, ls='--', color='r', linewidth=2, label=r'GTEMPO')

ax.scatter(betas2, data2, linewidth=0.5, marker='+', color='g', ls='--', label=r'CTMQC')

ax.set_title(r'$G({\tau})$,step=%s'%(step), fontsize=fontsize)

plt.tight_layout(pad=0.5)

# # plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


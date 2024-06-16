import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math


rc('text', usetex=True)


def read_rtempo(beta, t, i):
	filename = 'beta%st%s/Gtau-%s.dat'%(beta, t, i)
	data = loadtxt(filename)
	return data[:, 1]


def read_itempo(beta, i):
	mpath = '../imag/beta%s/'%(beta)
	filename = mpath + 'Gtau-%s.dat'%(i)
	data = loadtxt(filename)
	return data

fontsize = 20
labelsize = 16
linewidth = 2.
markersize = 10

colors = ['b', 'c', 'orange', 'g', 'y', 'r', 'k']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'r'


beta = 1
t = 10
step = 2

data1 = read_rtempo(beta, t, step)
betas1 = linspace(0, beta, num=len(data1))
data2 = read_itempo(beta, step)
betas2 = linspace(0, beta, num=len(data2))

fig, ax = plt.subplots(1,1, figsize=(8, 7))


ax.plot(betas1, data1, ls='--', color='r', linewidth=linewidth, label=r'rGTEMPO')

ax.plot(betas2, data2, ls='--', color='k', linewidth=linewidth, label=r'iGTEMPO')

ax.set_title(r'$G({\tau})$,step=%s'%(step), fontsize=fontsize)

ax.legend(fontsize=fontsize)

plt.tight_layout(pad=0.5)

# # plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


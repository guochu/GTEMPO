import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math


rc('text', usetex=True)


def read_tempo(beta, t, i):
	filename = 'beta%st%s/Gt-%s.dat'%(beta, t, i)
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

step = 1




fig, ax = plt.subplots(1,1, figsize=(8, 7))

t = 5
data1 = read_tempo(beta, t, step)
ax.plot(data1[:,0], data1[:, 2], ls='--', color='r', linewidth=2, label=r'$t=%s$'%(t))

t = 10
data1 = read_tempo(beta, t, step)
ax.plot(data1[:,0], data1[:, 2], ls='--', color='k', linewidth=2, label=r'$t=%s$'%(t))

ax.set_xlim(0, 20)

ax.set_title(r'$G({\tau})$,step=%s'%(step), fontsize=fontsize)

ax.legend(fontsize=12)

plt.tight_layout(pad=0.5)

# # plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


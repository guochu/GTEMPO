import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math


rc('text', usetex=True)


def read_tempo(beta, t, i):
	filename = 'beta%st%s/Aw-%s.dat'%(beta, t, i)
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

data1 = read_tempo(beta, t, step)


fig, ax = plt.subplots(1,1, figsize=(8, 7))


ax.plot(data1[:,0], data1[:, 1], ls='--', color='r', linewidth=2, label=r'GTEMPO')

ax.set_xlim(-3, 3)

plt.tight_layout(pad=0.5)

# # plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace
import math


rc('text', usetex=True)

def read_data(filename):
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], -asarray(data['gf'])

def read_mc_data(filename):
	data = loadtxt(filename)
	return data[:, 0], data[:, 1]


fontsize = 18
labelsize = 14
linewidth = 1.
markersize = 4

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

U = 1.
mu = U/2
beta = 60.
N = 600
chi = 512

data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s.json'%(beta, U, mu, N, chi)
ts, gf = read_data(data_path)

ax.plot(ts, gf, linewidth=linewidth, color='r', ls='-.', label=r'$\chi=%s, \beta=%s, N=%s$'%(chi, beta, N))




mc_data_path = '/Users/guochu/Documents/Since2018/GTEMPO/gtempo/chen/im/1orb/beta60/G-7.dat'
ts, gf = read_mc_data(mc_data_path)
ts = linspace(0, beta, num=len(ts))
ax.scatter(ts, gf, linewidth=0.5, marker='+', color='g', ls='--', label=r'CT-MQC')

ax.legend(fontsize=labelsize)

ax.set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax.set_xlabel(r'$\tau$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')


plt.tight_layout(pad=0.5)

plt.show()

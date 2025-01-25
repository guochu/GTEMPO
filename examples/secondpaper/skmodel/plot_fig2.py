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
	# print(data)
	return data[:, 0], data[:, 1]


fontsize = 18
labelsize = 14
linewidth = 1.
markersize = 4

fig, ax = plt.subplots(1, 2, figsize=(8, 4))

U = 2.
J = 0.5
mu = (3*U - 5*J)/2
# mu = U / 2
beta = 10.
N = 100


data_path = 'result/anderson_tempo1_norb2_beta%s_U%s_J%s_mu%s_N%s_b.json'%(beta, U, J, mu, N)
ts, gf = read_data(data_path)


ax[0].plot(ts, gf, linewidth=linewidth, color='r', ls='-.', markerfacecolor='none', label=r'GTEMPO')


mc_data_path = '/Users/guochu/Documents/Since2018/GTEMPO/gtempo/chen/im/2orbs/half/beta10/G-7.dat'
ts, gf = read_mc_data(mc_data_path)
ts = linspace(0, beta, num=len(ts))
ax[0].scatter(ts, gf, linewidth=0.5, color='g', marker='+', ls='--', label=r'CT-MQC')


ax[0].legend(fontsize=labelsize)

ax[0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].tick_params(axis='y', which='both')

mu = 2.5*U - 5*J
# mu = U / 2
beta = 2.
N = 20


data_path = 'result/anderson_tempo1_norb3_beta%s_U%s_J%s_mu%s_N%s_b.json'%(beta, U, J, mu, N)
ts, gf = read_data(data_path)


ax[1].plot(ts, gf, linewidth=linewidth, color='r', ls='-.', markerfacecolor='none', label=r'GTEMPO')


mc_data_path = '/Users/guochu/Documents/Since2018/GTEMPO/gtempo/chen/im/3orbs/half/beta2/G-7.dat'
ts, gf = read_mc_data(mc_data_path)
ts = linspace(0, beta, num=len(ts))
ax[1].scatter(ts, gf, linewidth=0.5, color='g', marker='+', ls='--', label=r'CT-MQC')


ax[1].legend(fontsize=labelsize)

ax[1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].tick_params(axis='y', which='both')

plt.tight_layout(pad=0.5)

# plt.savefig('multiorbital2.pdf', bbox_inches='tight')


plt.show()

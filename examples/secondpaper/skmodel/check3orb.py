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
	print(data)
	return data[:, 0], data[:, 1]


fontsize = 18
labelsize = 14
linewidth = 1.
markersize = 4

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

U = 2.
J = 0.5
mu = 2.5*U - 5*J
# mu = U / 2
beta = 2.
N = 20


data_path = 'result/anderson_tempo1_norb3_beta%s_U%s_J%s_mu%s_N%s_b.json'%(beta, U, J, mu, N)
ts, gf = read_data(data_path)

ax.plot(ts, gf, linewidth=linewidth, color='r', marker='o', ls='-.', markerfacecolor='none', label=r'GTEMPO')


chi = 60
chi2 = 240
chi3 = 1000
data_path = 'result/anderson_tempo1_norb3_beta%s_U%s_J%s_mu%s_N%s_chi%s_chi2%s_chi3%s.json'%(beta, U, J, mu, N, chi, chi2, chi3)
ts, gf = read_data(data_path)

ax.plot(ts, gf, linewidth=linewidth, color='b', marker='^', ls='-.', markerfacecolor='none', label=r'GTEMPO')


chi = 60
chi2 = 300
chi3 = 1000
data_path = 'result/anderson_tempo1_norb3_beta%s_U%s_J%s_mu%s_N%s_chi%s_chi2%s_chi3%s.json'%(beta, U, J, mu, N, chi, chi2, chi3)
ts, gf = read_data(data_path)

ax.plot(ts, gf, linewidth=linewidth, color='y', marker='+', ls='-.', markerfacecolor='none', label=r'GTEMPO')


mc_data_path = '/Users/guochu/Documents/Since2018/GTEMPO/gtempo/chen/im/3orbs/half/beta2/G-7.dat'
ts, gf = read_mc_data(mc_data_path)
ts = linspace(0, beta, num=len(ts))
ax.scatter(ts, gf, linewidth=0.5, color='g', marker='+', ls='--', label=r'CT-MQC')


ax.legend(fontsize=labelsize)

ax.set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax.set_xlabel(r'$\tau$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')


plt.tight_layout(pad=0.5)


# plt.savefig('multiorbital1.pdf', bbox_inches='tight')

plt.show()

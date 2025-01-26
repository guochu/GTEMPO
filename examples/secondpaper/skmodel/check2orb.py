import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
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

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)


fontsize = 18
labelsize = 14
linewidth = 1.
markersize = 4

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

U = 2.
J = 0.5
mu = (3*U - 5*J)/2
# mu = U / 2
beta = 10.
N = 100


data_path = 'result/anderson_tempo1_norb2_beta%s_U%s_J%s_mu%s_N%s_b.json'%(beta, U, J, mu, N)
ts, gf = read_data(data_path)

ax.plot(ts, gf, linewidth=linewidth, color='r', marker='o', ls='-.', markerfacecolor='none', label=r'Zipup')


chi = 80
chi2 = 320
chi3 = 1000
data_path = 'result/anderson_tempo1_norb2_beta%s_U%s_J%s_mu%s_N%s_chi%s_chi2%s_chi3%s.json'%(beta, U, J, mu, N, chi, chi2, chi3)
ts2, gf2 = read_data(data_path)

print('mse error ', mse_error(gf, gf2))

ax.plot(ts2, gf2, linewidth=linewidth, color='b', marker='^', ls='-.', markerfacecolor='none', label=r'PartialIntegration')



mc_data_path = '/Users/guochu/Documents/Since2018/GTEMPO/gtempo/chen/im/2orbs/half/beta%s/G-7.dat'%(round(beta))
ts, gf = read_mc_data(mc_data_path)
ts = linspace(0, beta, num=len(ts))
ax.scatter(ts, gf, linewidth=0.5, color='g', marker='+', ls='--', label=r'CTMQC')


ax.legend(fontsize=labelsize)

ax.set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax.set_xlabel(r'$\tau$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')


plt.tight_layout(pad=0.5)


# plt.savefig('multiorbital1.pdf', bbox_inches='tight')

plt.show()

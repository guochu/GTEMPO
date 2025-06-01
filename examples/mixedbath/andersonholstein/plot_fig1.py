import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math


rc('text', usetex=True)

def parse_complex_array(data):
	re = [item['re'] for item in data]
	im = [item['im'] for item in data]
	return asarray(re) + 1j * asarray(im)


def read_imag_tempo(beta, N, mu, d, alpha, chi=80):
	dt = beta / N
	filename = 'result/andersonholstein_imaggtempo_beta%s_dtau%s_d%s_alpha%s_mu%s_chi%s.json'%(beta, dt, d, alpha, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = data['gtau']
	gnn = data['nn']
	gnn2 = data['nn2']
	ts = asarray(data['taus'])
	return ts, gt, gnn, gnn2



def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 14
labelsize = 12
linewidth = 1.5
markersize = 4

colors = ['b', 'g', 'c', 'r', 'y']
markers = ['o', '^', '+', 'x', '*']

fig, ax = plt.subplots(1,2, figsize=(8,3.5))


# chi = 100

mu = 0.

beta = 5
N = 50
d = 3
# alpha = 0.5
alpha = 1


chis = [60, 80, 100, 120]

for i, chi in enumerate(chis):
	taus, gtau, gnn, gnn2 = read_imag_tempo(beta, N, mu, d, alpha, chi)
	ax[0].plot(taus, gtau, ls='--', color=colors[i],  markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))
	ax[1].plot(taus[:-1], gnn, ls='--', color=colors[i], markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))


ax[0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].legend(loc = 'center', fontsize=12)

ax[1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1].set_ylabel(r'$\langle \hat{n}(\tau)\hat{n}\rangle$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].locator_params(axis='both', nbins=6)
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

# plt.savefig('fig2.pdf', bbox_inches='tight')

plt.show()

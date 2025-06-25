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


def read_imag_tempo(beta, N, U, mu, d, alpha, chi=80):
	dt = beta / N
	filename = 'result/andersonholstein_int_imaggtempo_beta%s_dtau%s_d%s_alpha%s_U%s_mu%s_chi%s.json'%(beta, dt, d, alpha, U, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = data['gtau']
	gnn = data['nn']
	gnn.append(gnn[0])
	ts = asarray(data['taus'])
	return ts-ts[0], gt, gnn



def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 20
labelsize = 16
linewidth = 2
markersize = 8

colors = ['b', 'g', 'c', 'r', 'y']
markers = ['o', '^', '+', 'x', '*']

fig, ax = plt.subplots(1,2, figsize=(8,4))

annotate_xy = (-0.15, 1.07)

# chi = 100

U = 1
mu = U/2


chi = 100	

beta = 5
Ns = [50, 100]
d = 1
alpha = 0.1

for (i, N) in enumerate(Ns):
	taus, gtau, gnn = read_imag_tempo(beta, N, U, mu, d, alpha, chi)
	print('chi=', chi, ' ', gtau[-1], ' ', gnn[0])
	ax[0].plot(taus, gtau, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\delta=%s$'%(beta/N))
	ax[1].plot(taus, gnn, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\delta=%s$'%(beta/N))


ax[0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].annotate(r'(a1)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0].legend( fontsize=12)

ax[1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].locator_params(axis='both', nbins=6)
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1].annotate(r'(b1)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[1].legend( fontsize=12)





plt.tight_layout(pad=0.5)

# plt.savefig('full_noint_imag.pdf', bbox_inches='tight')

plt.show()

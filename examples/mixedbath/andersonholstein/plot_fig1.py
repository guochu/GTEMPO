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
linewidth = 3
markersize = 8


linewidth_s = 2
markersize_s = 4
fontsize_s = 16
labelsize_s = 14

colors = ['b', 'g', 'c', 'r', 'y']
markers = ['o', '^', '+', 'x', '*']

fig, ax = plt.subplots(2,2, figsize=(8,7))

annotate_xy = (-0.15, 1.07)

# chi = 100

mu = 0

beta = 1
N_max = 40
d = 3
alpha = 1

chi = 100

taus, gtau, gnn = read_imag_tempo(beta, N_max, mu, d, alpha, chi)

ax[0,0].plot(taus, gtau, ls='--', color='k', markersize=markersize, markerfacecolor='none', linewidth=linewidth)

ax[1,0].plot(taus, gnn, ls='--', color='b', markersize=markersize, markerfacecolor='none', linewidth=linewidth)


ax[0,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,0].set_title(r'$\beta=%s$'%(beta), fontsize=labelsize)

ax[1,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(c)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


Ns = [5, 10, 20]

dtaus = [beta / N for N in Ns]


gtau_errs = []
nn_errs = []


for i, N in enumerate(Ns):
	taus1, gtau1, gnn1 = read_imag_tempo(beta, N, mu, d, alpha, chi)

	step = N_max // N
	gtau_scaled = gtau[0:step:len(gtau)]
	gnn_scaled = gnn[0:step:len(gnn)]

	gtau_errs.append(mse_error(gtau_scaled, gtau1[:len(gtau_scaled)]))
	nn_errs.append(mse_error(gnn_scaled, gnn1[:len(gnn_scaled)]))


ax1 = ax[0,0].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.semilogy(dtaus, gtau_errs, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)


ax1 = ax[1,0].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.semilogy(dtaus, nn_errs, ls='--', color='b', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)


# beta = 10
chi = 200
beta = 10
N_max = 400


taus, gtau, gnn = read_imag_tempo(beta, N_max, mu, d, alpha, chi)

ax[0,1].plot(taus, gtau, ls='--', color='k', markersize=markersize, markerfacecolor='none', linewidth=linewidth)

ax[1,1].plot(taus, gnn, ls='--', color='b', markersize=markersize, markerfacecolor='none', linewidth=linewidth)


ax[0,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,1].set_title(r'$\beta=%s$'%(beta), fontsize=labelsize)

ax[1,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(d)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


Ns = [50, 100, 200]

dtaus = [beta / N for N in Ns]


gtau_errs = []
nn_errs = []


for i, N in enumerate(Ns):
	taus1, gtau1, gnn1 = read_imag_tempo(beta, N, mu, d, alpha, chi)

	step = N_max // N
	gtau_scaled = gtau[0:step:len(gtau)]
	gnn_scaled = gnn[0:step:len(gnn)]

	gtau_errs.append(mse_error(gtau_scaled, gtau1[:len(gtau_scaled)]))
	nn_errs.append(mse_error(gnn_scaled, gnn1[:len(gnn_scaled)]))


ax1 = ax[0,1].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.semilogy(dtaus, gtau_errs, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)


ax1 = ax[1,1].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.semilogy(dtaus, nn_errs, ls='--', color='b', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)

plt.tight_layout(pad=0.5)

plt.savefig('full_imag_noint.pdf', bbox_inches='tight')

plt.show()

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

fig, ax = plt.subplots(2,4, figsize=(12,6))

annotate_xy = (-0.15, 1.07)

# chi = 100

mu = 0

beta = 1
N = 20
d = 3
alpha = 1


chis = [20, 40, 60, 100]


for i, chi in enumerate(chis):
	taus, gtau, gnn = read_imag_tempo(beta, N, mu, d, alpha, chi)
	print('chi=', chi, ' ', gtau[-1], ' ', gnn[0])
	ax[0,0].plot(taus, gtau, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))
	ax[0,1].plot(taus[:-1], gnn, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))


ax[0,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a1)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend( fontsize=12)

ax[0,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b1)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,1].legend( fontsize=12)

Ns = [5, 10, 20, 40, 80]


chi = 60
for i, N in enumerate(Ns):
	taus, gtau, gnn = read_imag_tempo(beta, N, mu, d, alpha, chi)
	ax[0,2].plot(taus, gtau, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\delta t=%s$'%(beta/N))
	ax[0,3].plot(taus[:-1], gnn, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\delta t=%s$'%(beta/N))

ax[0,2].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,2].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,2].locator_params(axis='both', nbins=6)
ax[0,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,2].annotate(r'(c1)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,2].legend(fontsize=12)

ax[0,3].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,3].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[0,3].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,3].locator_params(axis='both', nbins=6)
ax[0,3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,3].annotate(r'(d1)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,3].legend(fontsize=12)


# beta = 10
chis = [120, 140, 160]
beta = 10
N = 100

for i, chi in enumerate(chis):
	taus, gtau, gnn = read_imag_tempo(beta, N, mu, d, alpha, chi)
	print('chi=', chi, ' ', gtau[-1], ' ', gnn[0])
	ax[1,0].plot(taus, gtau, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))
	ax[1,1].plot(taus[:-1], gnn, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))


ax[1,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(a2)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[1,0].legend( fontsize=12)

ax[1,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(b2)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[1,1].legend( fontsize=12)

Ns = [50, 100, 200]

chi = 200
for i, N in enumerate(Ns):
	taus, gtau, gnn = read_imag_tempo(beta, N, mu, d, alpha, chi)
	ax[1,2].plot(taus, gtau, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\delta t=%s$'%(beta/N))
	ax[1,3].plot(taus[:-1], gnn, ls='--', color=colors[i], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'$\delta t=%s$'%(beta/N))

ax[1,2].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,2].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,2].locator_params(axis='both', nbins=6)
ax[1,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,2].annotate(r'(c2)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[1,2].legend(fontsize=12)

ax[1,3].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,3].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1,3].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,3].locator_params(axis='both', nbins=6)
ax[1,3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,3].annotate(r'(d2)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[1,3].legend(fontsize=12)


plt.tight_layout(pad=0.5)

# plt.savefig('full_noint_imag.pdf', bbox_inches='tight')

plt.show()

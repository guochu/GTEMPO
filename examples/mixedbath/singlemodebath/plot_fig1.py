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


def read_imag_tempo(beta, Ntau, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dtau = beta / Ntau
	filename = 'result/noninteracting_imaggtempo_beta%s_dtau%s_omega0%s_alpha0%s_omega1%s_alpha1%s_mu%s_chi%s.json'%(beta, dtau, omega0, alpha0, omega1, alpha1, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gf = asarray(data['gf'])
	ts = asarray(data['taus'])
	nn = data['nn']
	nn.append(nn[0])
	return ts, gf, nn


def read_imag_ed(beta, Ntau, mu, omega0, alpha0, omega1, alpha1):
	dtau = beta / Ntau
	filename = 'result/noninteracting_eq_ED_imag_beta%s_mu%s_dtau%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, dtau, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray(data['gt'])
	ts = asarray(data['taus'])
	nn = asarray(data['nn'])
	return ts, gt, nn


def read_imag_tempo_int(beta, Ntau, U, J, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dtau = beta / Ntau
	filename = 'result/interacting_imaggtempo_beta%s_dtau%s_omega0%s_alpha0%s_omega1%s_alpha1%s_U%s_J%s_mu%s_chi%s.json'%(beta, dtau, omega0, alpha0, omega1, alpha1, U, J, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gf = asarray(data['gf'])
	ts = asarray(data['taus'])
	nn = data['nn']
	nn.append(nn[0])
	return ts, gf, nn


def read_imag_ed_int(beta, Ntau, U, J, mu, omega0, alpha0, omega1, alpha1):
	dtau = beta / Ntau
	filename = 'result/interacting_eq_ED_imag_beta%s_U%s_J%s_mu%s_dtau%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, U, J, mu, dtau, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray(data['gt'])
	ts = asarray(data['taus'])
	nn = asarray(data['nn'])
	return ts, gt, nn

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 20
labelsize = 16
linewidth1 = 1.5
linewidth2 = 3
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2,2, figsize=(8,6.5))


omega0 = 1
alpha0 = 0.5
omega1 = 1
alpha1 = 1
chi = 50

mu = 0

beta = 5
Ntau = 200

# real time data
ts, gt, nn = read_imag_ed(beta, Ntau, mu, omega0, alpha0, omega1, alpha1)

ax[0,0].plot(ts, gt.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')
ax[0,1].plot(ts, nn.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')

ts2, gt2, nn2 = read_imag_tempo(beta, Ntau, mu, omega0, alpha0, omega1, alpha1, chi)

ax[0,0].plot(ts2, gt2.real, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none', label=r'$\chi=%s$'%(chi))
ax[0,1].plot(ts2, nn2, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none', label=r'$\chi=%s$'%(chi))


ax[0,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[0,0].legend(loc='lower right', fontsize=12)

ax[0,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


U = 2
J = 1
mu = U/2


ts, gt, nn = read_imag_ed_int(beta, Ntau, U, J, mu, omega0, alpha0, omega1, alpha1)
ts2, gt2, nn2 = read_imag_tempo_int(beta, Ntau, U, J, mu, omega0, alpha0, omega1, alpha1, chi)

ax[1,0].plot(ts, gt.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')
ax[1,1].plot(ts, nn.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')

ax[1,0].plot(ts2, gt2.real, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none', label=r'$\chi=%s$'%(chi))
ax[1,1].plot(ts2, nn2, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none', label=r'$\chi=%s$'%(chi))

ax[1,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(c)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[0,0].legend(loc='lower right', fontsize=12)

ax[1,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(d)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

# plt.savefig('independentbosons1.pdf', bbox_inches='tight')

plt.show()

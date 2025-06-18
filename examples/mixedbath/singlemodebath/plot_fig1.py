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
	return ts, gf


def read_imag_ed(beta, Ntau, mu, omega0, alpha0, omega1, alpha1):
	dtau = beta / Ntau
	filename = 'result/noninteracting_eq_ED_imag_beta%s_mu%s_dtau%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, dtau, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray(data['gt'])
	ts = asarray(data['taus'])
	return ts, gt


def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 20
labelsize = 16
linewidth = 2.5
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(1,2, figsize=(8,3.5))


omega0 = 1
alpha0 = 0.5
omega1 = 1
alpha1 = 1
chi = 20

mu = 0.

beta = 0.5
Ntau = 10

# real time data
ts, gt = read_imag_ed(beta, Ntau, mu, omega0, alpha0, omega1, alpha1)

ax[0].plot(ts, gt.real, ls='--', color='k', linewidth=linewidth, label=r'ED')

ts2, gt2 = read_imag_tempo(beta, Ntau, mu, omega0, alpha0, omega1, alpha1, chi)


ax[0].plot(ts2, gt2.real, color='k', ls='none', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none', label=r'$\chi=%s$'%(chi))


ax[0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].legend(loc='lower right', fontsize=12)


errs = []
chis = [4,8,12,16,20]
for i, chi in enumerate(chis):
	ts2, gt2 = read_imag_tempo(beta, Ntau, mu, omega0, alpha0, omega1, alpha1, chi)
	errs.append(mse_error(gt, gt2))


ax[1].plot(chis, errs, ls='--', color='c', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none')

ax[1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].locator_params(axis='both', nbins=6)
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[1].legend(loc='lower right', fontsize=12)


plt.tight_layout(pad=0.5)

# plt.savefig('independentbosons1.pdf', bbox_inches='tight')

plt.show()

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

linewidth_s = 2
markersize_s = 4
fontsize_s = 16
labelsize_s = 14


colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2,2, figsize=(8,7))


omega0 = 1
alpha0 = 0.5
omega1 = 1
alpha1 = 1
chi = 100

annotate_xy = (-0.15, 1.05)

mu = 0

beta = 5
Ntau = 200

Ntaus = [25, 50, 100, 200]

# real time data
ts, gt, nn = read_imag_ed(beta, Ntau, mu, omega0, alpha0, omega1, alpha1)

ax[0,0].plot(ts, gt.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')
ax[1,0].plot(ts, nn.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')


ts2, gt2, nn2 = read_imag_tempo(beta, Ntau, mu, omega0, alpha0, omega1, alpha1, chi)
ax[0,0].plot(ts2, gt2.real, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none')
ax[1,0].plot(ts2, nn2, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none')


ax[0,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,0].set_title(r'Single flavor', fontsize=labelsize)

ax[1,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(c)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


gt_errs = []
nn_errs = []
dtaus = [beta/Nt for Nt in Ntaus]

for i, Nt in enumerate(Ntaus):
	dtau = beta / Nt
	ts, gt, nn = read_imag_ed(beta, Nt, mu, omega0, alpha0, omega1, alpha1)
	ts2, gt2, nn2 = read_imag_tempo(beta, Nt, mu, omega0, alpha0, omega1, alpha1, chi)
	gt_errs.append(mse_error(gt, gt2))
	nn_errs.append(mse_error(nn, nn2))

ax1 = ax[0,0].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.plot(dtaus, gt_errs, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')


ax1.set_ylabel(r'$\mathcal{E}[G(\tau)]$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax1 = ax[1,0].inset_axes([0.3, 0.4, 0.5, 0.5])

ax1.plot(dtaus, nn_errs, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')


ax1.set_ylabel(r'$\mathcal{E}[X(\tau)]$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


# two flavor
U = 2
J = 1
mu = U/2


ts, gt, nn = read_imag_ed_int(beta, Ntau, U, J, mu, omega0, alpha0, omega1, alpha1)

ax[0,1].plot(ts, gt.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')
ax[1,1].plot(ts, nn.real, ls='-', color='k', linewidth=linewidth1, label=r'ED')


ts2, gt2, nn2 = read_imag_tempo_int(beta, Ntau, U, J, mu, omega0, alpha0, omega1, alpha1, chi)
ax[0,1].plot(ts2, gt2.real, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none', label=r'$\delta\tau=%s$'%(dtau))
ax[1,1].plot(ts2, nn2, color='k', ls='--', linewidth=linewidth2, markersize=markersize, markerfacecolor='none', label=r'$\delta\tau=%s$'%(dtau))


ax[0,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
ax[0,1].set_title(r'Two flavor', fontsize=labelsize)

ax[1,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$X(\tau)$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(d)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

gt_errs = []
nn_errs = []
# dtaus = [beta/Nt for Nt in Ntaus]

for i, Nt in enumerate(Ntaus):
	dtau = beta / Nt
	ts, gt, nn = read_imag_ed_int(beta, Nt, U, J, mu, omega0, alpha0, omega1, alpha1)
	ts2, gt2, nn2 = read_imag_tempo_int(beta, Nt, U, J, mu, omega0, alpha0, omega1, alpha1, chi)
	gt_errs.append(mse_error(gt, gt2))
	nn_errs.append(mse_error(nn, nn2))

ax1 = ax[0,1].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.plot(dtaus, gt_errs, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')


ax1.set_ylabel(r'$\mathcal{E}[G(\tau)]$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax1 = ax[1,1].inset_axes([0.3, 0.4, 0.5, 0.5])

ax1.plot(dtaus, nn_errs, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')


ax1.set_ylabel(r'$\mathcal{E}[X(\tau)]$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\delta\tau$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


plt.tight_layout(pad=0.5)

plt.savefig('toy_imag.pdf', bbox_inches='tight')

plt.show()

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

def read_real_tempo(beta, t, N, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dt = t / N
	filename = 'result/noninteracting_realgtempo_beta%s_t%s_dt%s_omega0%s_alpha0%s_omega1%s_alpha1%s_mu%s_chi%s.json'%(beta, t, dt, omega0, alpha0, omega1, alpha1, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt, nn

def read_mixed_tempo(beta, Ntau, t, N, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dt = t / N
	dtau = beta / Ntau
	filename = 'result/noninteracting_mixedgtempo_beta%s_dtau%s_t%s_dt%s_omega0%s_alpha0%s_omega1%s_alpha1%s_mu%s_chi%s.json'%(beta, dtau, t, dt, omega0, alpha0, omega1, alpha1, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt, nn

def read_neq_ed(beta, t, N, mu, omega0, alpha0, omega1, alpha1):
	filename = 'result/noninteracting_neq_ED_real_beta%s_mu%s_t%s_N%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, t, N, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	return data['ts'], gt, lt, nn

def read_eq_ed(beta, t, N, mu, omega0, alpha0, omega1, alpha1):
	filename = 'result/noninteracting_eq_ED_real_beta%s_mu%s_t%s_N%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, t, N, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	return data['ts'], gt, lt, nn


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

fig, ax = plt.subplots(3,3, figsize=(11,9))


omega0 = 1
alpha0 = 0.5
omega1 = 1
alpha1 = 1
chi = 600

mu = 0.

t = 5
Nt = 200
beta = 5
Ntau = 100

# mixed time data
ts, gt, lt, nn = read_eq_ed(beta, t, Nt, mu, omega0, alpha0, omega1, alpha1)
ts2, gt2, lt2, nn2 = read_mixed_tempo(beta, Ntau, t, Nt, mu, omega0, alpha0, omega1, alpha1, chi)

annotate_xy = (-0.15, 1.07)

ax1color = 'purple'

ax[0,0].plot(ts, gt.real, ls='-', color=ax1color, linewidth=linewidth1, label=r'ED')
ax[0,0].plot(ts2, gt2.real, ls='--', color=ax1color, linewidth=linewidth2, label=r'GTEMPO')


ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize, color=ax1color)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize, colors=ax1color)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
# ax[0,0].legend(loc='upper right', fontsize=12)


ax2color = 'c'

ax2 = ax[0,0].twinx()

ax2.plot(ts, gt.imag, ls='-', color=ax2color, linewidth=linewidth1, label=r'ED')
ax2.plot(ts2, gt2.imag, ls='--', color=ax2color, linewidth=linewidth2, label=r'GTEMPO')

ax2.set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))



ax[0,1].plot(ts, lt.real, ls='-', color=ax1color, linewidth=linewidth1, label=r'ED')
ax[0,1].plot(ts2, lt2.real, ls='--', color=ax1color, linewidth=linewidth2, label=r'GTEMPO')

ax[0,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,1].set_ylabel(r'${\rm Re}[G^{<}(t)]$', fontsize=fontsize, color=ax1color)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize, colors=ax1color)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


ax2 = ax[0,1].twinx()

ax2.plot(ts, lt.imag, ls='-', color=ax2color, linewidth=linewidth1, label=r'ED')
ax2.plot(ts2, lt2.imag, ls='--', color=ax2color, linewidth=linewidth2, label=r'GTEMPO')


ax2.set_ylabel(r'${\rm Im}[G^{<}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax3color = 'b'

ax[0,2].plot(ts[2:], nn.real[2:], ls='-', color=ax3color, linewidth=linewidth1, label=r'ED')
ax[0,2].plot(ts2[1:-1], nn2.real, ls='--', color=ax3color, linewidth=linewidth2, markersize=markersize, markerfacecolor='none', label=r'GTEMPO, $\chi=%s$'%(chi))

ax[0,2].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,2].set_ylabel(r'$X(t)$', fontsize=fontsize)
ax[0,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,2].locator_params(axis='both', nbins=6)
ax[0,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,2].annotate(r'(c)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


chis = [100, 200, 300, 400, 500, 600, 700]
nn_errs = []
gt_errs = []
lt_errs = []

for chi in chis:
	ts2, gt2, lt2, nn2 = read_mixed_tempo(beta, Ntau, t, Nt, mu, omega0, alpha0, omega1, alpha1, chi)
	gt_errs.append(mse_error(gt, gt2))
	lt_errs.append(mse_error(lt, lt2))
	nn_errs.append(mse_error(nn.real[1:-1], nn2.real))

ax[1,0].plot(chis, gt_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[1,0].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(d)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[1,1].plot(chis, lt_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[1,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(e)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[1,2].plot(chis, nn_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[1,2].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,2].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,2].locator_params(axis='both', nbins=6)
ax[1,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,2].annotate(r'(f)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)



Nts = [50, 100, 200]
dts = [t / Nt for Nt in Nts]
chi = 600
gt_errs = []
lt_errs = []
nn_errs = []

for Nt in Nts:
	ts, gt, lt, nn = read_eq_ed(beta, t, Nt, mu, omega0, alpha0, omega1, alpha1)
	ts2, gt2, lt2, nn2 = read_mixed_tempo(beta, Ntau, t, Nt, mu, omega0, alpha0, omega1, alpha1, chi)
	gt_errs.append(mse_error(gt, gt2))
	lt_errs.append(mse_error(lt, lt2))
	nn_errs.append(mse_error(nn.real[1:-1], nn2.real))


ax[2,0].plot(dts, gt_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[2,0].set_xlabel(r'$\delta t$', fontsize=fontsize)
ax[2,0].set_ylabel(r'$\mathcal{E}[G^{>}(t)]$', fontsize=fontsize)
ax[2,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,0].locator_params(axis='both', nbins=6)
ax[2,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,0].annotate(r'(g)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[2,1].plot(dts, lt_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[2,1].set_xlabel(r'$\delta t$', fontsize=fontsize)
ax[2,1].set_ylabel(r'$\mathcal{E}[G^{<}(t)]$', fontsize=fontsize)
ax[2,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,1].locator_params(axis='both', nbins=6)
ax[2,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,1].annotate(r'(h)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[2,2].plot(dts, nn_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[2,2].set_xlabel(r'$\delta t$', fontsize=fontsize)
ax[2,2].set_ylabel(r'$\mathcal{E}[X(t)]$', fontsize=fontsize)
ax[2,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,2].locator_params(axis='both', nbins=6)
ax[2,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,2].annotate(r'(i)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

plt.savefig('toy_mixed.pdf', bbox_inches='tight')

plt.show()

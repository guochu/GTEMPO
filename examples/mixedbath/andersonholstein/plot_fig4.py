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


def read_real_tempo(beta, t, N, U, mu, d, alpha, chi=80):
	dt = t / N
	filename = 'result/andersonholstein_int_realgtempo_beta%s_t%s_dt%s_d%s_alpha%s_U%s_mu%s_chi%s.json'%(beta, t, dt, d, alpha, U, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gnn = parse_complex_array(data['nn'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt, gnn



def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 20
labelsize = 16
linewidth = 3
markersize = 10

colors = ['b', 'g', 'c', 'r', 'y']
markers = ['o', '^', '+', 'x', '*']

fig, ax = plt.subplots(3,3, figsize=(11, 9))

annotate_xy = (-0.15, 1.07)

ax1color = 'purple'
ax2color = 'c'
ax3color = 'b'


# chi = 100

U = 1
mu = U / 2

beta = 10
t = 1

chi_max = 400
N_max = 40

d = 3
alpha = 1

ts, gt, lt, gnn = read_real_tempo(beta, t, N_max, U, mu, d, alpha, chi_max)

ax[0,0].plot(ts, gt.real, ls='--', color=ax1color, markerfacecolor='none', linewidth=linewidth)
ax[1,0].plot(ts, lt.real, ls='--', color=ax1color, markerfacecolor='none', linewidth=linewidth)
ax[2,0].plot(ts[:-1], gnn.real, ls='--', color=ax3color, markerfacecolor='none', linewidth=linewidth)


ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize, color=ax1color)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].tick_params(axis='y', colors=ax1color)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
# ax[0,0].legend(loc = 'center', fontsize=12)

ax[1,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,0].set_ylabel(r'${\rm Re}[G^{<}(t)]$', fontsize=fontsize, color=ax1color)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].tick_params(axis='y', colors=ax1color)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(d)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[2,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[2,0].set_ylabel(r'$X(t)$', fontsize=fontsize)
ax[2,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,0].locator_params(axis='both', nbins=6)
ax[2,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,0].annotate(r'(g)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


ax2 = ax[0,0].twinx()
ax2.plot(ts, gt.imag, ls='--', color=ax2color, linewidth=linewidth, label=r'ED')


ax2.set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


ax2 = ax[1,0].twinx()
ax2.plot(ts, lt.imag, ls='--', color=ax2color, linewidth=linewidth, label=r'ED')


ax2.set_ylabel(r'${\rm Im}[G^{<}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


chis = [100, 200, 300]

gt_errs = []
lt_errs = []
nn_errs = []

for i, chi in enumerate(chis):
	ts1, gt1, lt1, gnn1 = read_real_tempo(beta, t, N_max, U, mu, d, alpha, chi)
	gt_errs.append(mse_error(gt, gt1))
	lt_errs.append(mse_error(lt, lt1))
	nn_errs.append(mse_error(gnn, gnn1))


ax[0,1].plot(chis, gt_errs, ls='--', color='k', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none')
ax[0,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$\mathcal{E}[G^{>}(t)]$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[1,1].plot(chis, lt_errs, ls='--', color='k', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none')
ax[1,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$\mathcal{E}[G^{<}(t)]$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(e)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[2,1].plot(chis, nn_errs, ls='--', color='k', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none')
ax[2,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[2,1].set_ylabel(r'$\mathcal{E}[X(t)]$', fontsize=fontsize)
ax[2,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,1].locator_params(axis='both', nbins=6)
ax[2,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,1].annotate(r'(h)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


Nts = [5, 10, 20]
dts = [t / Nt for Nt in Nts]
# chi = 200

gt_errs = []
lt_errs = []
nn_errs = []


for i, Nt in enumerate(Nts):
	ts1, gt1, lt1, gnn1 = read_real_tempo(beta, t, Nt, U, mu, d, alpha, chi_max)
	step = N_max // Nt
	gt_scaled = gt[0:step:len(gt)]
	lt_scaled = lt[0:step:len(lt)]
	gnn_scaled = gnn[0:step:len(gnn)]
	gt_errs.append(mse_error(gt_scaled, gt1[:len(gt_scaled)]))
	lt_errs.append(mse_error(lt_scaled, lt1[:len(lt_scaled)]))
	nn_errs.append(mse_error(gnn_scaled, gnn1[:len(gnn_scaled)]))


ax[0,2].plot(dts, gt_errs, ls='--', color='k', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none')
ax[0,2].set_xlabel(r'$\delta t$', fontsize=fontsize)
ax[0,2].set_ylabel(r'$\mathcal{E}[G^{>}(t)]$', fontsize=fontsize)
ax[0,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,2].locator_params(axis='both', nbins=6)
ax[0,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,2].annotate(r'(c)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[1,2].plot(dts, lt_errs, ls='--', color='k', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none')
ax[1,2].set_xlabel(r'$\delta t$', fontsize=fontsize)
ax[1,2].set_ylabel(r'$\mathcal{E}[G^{<}(t)]$', fontsize=fontsize)
ax[1,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,2].locator_params(axis='both', nbins=6)
ax[1,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,2].annotate(r'(f)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


ax[2,2].plot(dts, nn_errs, ls='--', color='k', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none')
ax[2,2].set_xlabel(r'$\delta t$', fontsize=fontsize)
ax[2,2].set_ylabel(r'$\mathcal{E}[X(t)]$', fontsize=fontsize)
ax[2,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,2].locator_params(axis='both', nbins=6)
ax[2,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,2].annotate(r'(i)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

plt.savefig('full_real_int.pdf', bbox_inches='tight')

plt.show()

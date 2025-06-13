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
	gnn2 = parse_complex_array(data['nn2'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt, gnn, gnn2



def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 20
labelsize = 16
linewidth = 2
markersize = 10

colors = ['b', 'g', 'c', 'r', 'y']
markers = ['o', '^', '+', 'x', '*']

fig, ax = plt.subplots(2,3, figsize=(10, 5.5))

annotate_xy = (-0.15, 1.07)

ax1color = 'purple'
ax2color = 'c'
ax3color = 'b'


# chi = 100

U = 1
mu = U / 2

beta = 10
t = 1
N = 20
d = 1
alpha = 0.1


chis = [100, 200, 300]

alphas = [0.2, 0.5, 0.7, 1]

for i, chi in enumerate(chis):
	ts, gt, lt, gnn, gnn2 = read_real_tempo(beta, t, N, U, mu, d, alpha, chi)
	ax[0,0].plot(ts, gt.real, ls='--', color=ax1color, alpha=alphas[i], markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))
	ax[0,1].plot(ts, lt.real, ls='--', color=ax1color, alpha=alphas[i], markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))
	ax[0,2].plot(ts[:-1], gnn.real, ls='--', color=ax3color, alpha=alphas[i], markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))


ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize, color=ax1color)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].tick_params(axis='y', colors=ax1color)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
# ax[0,0].legend(loc = 'center', fontsize=12)

ax[0,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,1].set_ylabel(r'${\rm Re}[G^{<}(t)]$', fontsize=fontsize, color=ax1color)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].tick_params(axis='y', colors=ax1color)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[0,2].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,2].set_ylabel(r'$X(t)$', fontsize=fontsize)
ax[0,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,2].locator_params(axis='both', nbins=6)
ax[0,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,2].annotate(r'(c)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


ax2 = ax[0,0].twinx()
for i, chi in enumerate(chis):
	ts, gt, lt, gnn, gnn2 = read_real_tempo(beta, t, N, U, mu, d, alpha, chi)
	ax2.plot(ts, gt.imag, ls='--', color=ax2color, alpha=alphas[i], linewidth=linewidth, label=r'ED')


ax2.set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


ax2 = ax[0,1].twinx()
for i, chi in enumerate(chis):
	ts, gt, lt, gnn, gnn2 = read_real_tempo(beta, t, N, U, mu, d, alpha, chi)
	ax2.plot(ts, gt.imag, ls='--', color=ax2color, alpha=alphas[i], linewidth=linewidth, label=r'ED')


ax2.set_ylabel(r'${\rm Im}[G^{<}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


Nts = [5, 10, 20]
dts = [t / Nt for Nt in Nts]
chi = 300


for i, Nt in enumerate(Nts):
	ts, gt, lt, gnn, gnn2 = read_real_tempo(beta, t, Nt, U, mu, d, alpha, chi)
	ax[1,0].plot(ts, gt.real, ls='--', color=ax1color, alpha=alphas[i], linewidth=linewidth, label=r'$\chi=%s$'%(chi))
	ax[1,1].plot(ts, lt.real, ls='--', color=ax1color, alpha=alphas[i], linewidth=linewidth, label=r'$\chi=%s$'%(chi))
	ax[1,2].plot(ts[:-1], gnn.real, ls='--', color=ax3color, alpha=alphas[i], markerfacecolor='none', linewidth=linewidth, label=r'$\chi=%s$'%(chi))


ax[1,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize, color=ax1color)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].tick_params(axis='y', colors=ax1color)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(d)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)
# ax[0,0].legend(loc = 'center', fontsize=12)

ax[1,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,1].set_ylabel(r'${\rm Re}[G^{<}(t)]$', fontsize=fontsize, color=ax1color)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].tick_params(axis='y', colors=ax1color)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(e)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


ax[1,2].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,2].set_ylabel(r'$X(t)$', fontsize=fontsize)
ax[1,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,2].locator_params(axis='both', nbins=6)
ax[1,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,2].annotate(r'(f)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


ax2 = ax[1,0].twinx()
for i, Nt in enumerate(Nts):
	ts, gt, lt, gnn, gnn2 = read_real_tempo(beta, t, Nt, U, mu, d, alpha, chi)
	ax2.plot(ts, gt.imag, ls='--', color=ax2color, alpha=alphas[i], linewidth=linewidth, label=r'ED')


ax2.set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


ax2 = ax[1,1].twinx()
for i, Nt in enumerate(Nts):
	ts, gt, lt, gnn, gnn2 = read_real_tempo(beta, t, Nt, U, mu, d, alpha, chi)
	ax2.plot(ts, gt.imag, ls='--', color=ax2color, alpha=alphas[i], linewidth=linewidth, label=r'ED')


ax2.set_ylabel(r'${\rm Im}[G^{<}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

plt.tight_layout(pad=0.5)

# plt.savefig('fig2.pdf', bbox_inches='tight')

plt.show()

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

def read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi=80):
	dtau = beta / Ntau
	filename = 'result/noninteracting_imaggtempo_beta%s_dtau%s_mu%s_d%s_alpha%s_chi%s.json'%(beta, dtau, mu, d, alpha, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gtau = asarray(data['gtau'])
	taus = asarray(data['taus'])
	return taus-taus[0], gtau

def read_noninteracting_mixed_tempo(beta, Ntau, t, N, mu, d, alpha, chi=80):
	dt = t / N
	dtau = beta / Ntau
	filename = 'result/noninteracting_mixedgtempo_beta%s_dtau%s_t%s_dt%s_mu%s_d%s_alpha%s_chi%s.json'%(beta, dtau, t, dt, mu, d, alpha, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gtau = parse_complex_array(data['gtau'])
	ts = asarray(data['ts'])
	taus = asarray(data['taus'])
	return ts-ts[0], taus-taus[0], gt, lt, gtau

def read_noninteracting_imag_analytic(beta, N, mu, d, alpha):
	filename = 'result/noninteracting_analytic_imag_beta%s_mu%s_N%s_d%s_alpha%s.json'%(beta, mu, N, d, alpha)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = data['gf']
	return data['ts'], gt

def read_noninteracting_real_analytic(beta, t, N, mu, d, alpha):
	filename = 'result/noninteracting_analytic_real_beta%s_mu%s_t%s_N%s_d%s_alpha%s.json'%(beta, mu, t, N, d, alpha)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	return data['ts'], gt, lt


def read_interacting_mixed_tempo(beta, Ntau, t, N, U, mu, d, alpha, chi=80):
	dt = t / N
	dtau = beta / Ntau
	filename = 'result/interacting_mixedgtempo_beta%s_dtau%s_t%s_dt%s_U%s_mu%s_d%s_alpha%s_chi%s.json'%(beta, dtau, t, dt, U, mu, d, alpha, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gtau = parse_complex_array(data['gtau'])
	ts = asarray(data['ts'])
	taus = asarray(data['taus'])
	return ts-ts[0], taus-taus[0], gt, lt, gtau

def read_interacting_imag_analytic(beta, N, U, mu, d, alpha):
	filename = 'result/interacting_analytic_imag_beta%s_U%s_mu%s_N%s_d%s_alpha%s.json'%(beta, U, mu, N, d, alpha)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = data['gf']
	return data['ts'], gt

def read_interacting_real_analytic(beta, t, N, U, mu, d, alpha):
	filename = 'result/interacting_analytic_real_beta%s_U%s_mu%s_t%s_N%s_d%s_alpha%s.json'%(beta, U, mu, t, N, d, alpha)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	return data['ts'], gt, lt

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 22
labelsize = 18
linewidth1 = 1.5
linewidth2 = 3
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2,2, figsize=(8,6.5))


chi = 140

mu = 0.5

t = 5
Nt = 100
beta = 5
Ntau = 50
d = 1
alpha = 1

# noninteracting case
# taus, gtau = read_noninteracting_imag_analytic(beta, Ntau, mu, d, alpha)
ts, gt, lt = read_noninteracting_real_analytic(beta, t, Nt, mu, d, alpha)
# gf = gt - lt

ts2, taus2, gt2, lt2, gtau2 = read_noninteracting_mixed_tempo(beta, Ntau, t, Nt, mu, d, alpha, chi)
# gf2 = gt2 - lt2

annotate_xy = (-0.15, 1.05)

ax1color = 'purple'

ax[0,0].plot(ts, gt.real, ls='-', color=ax1color, linewidth=linewidth1, label=r'Analytic')
ax[0,0].plot(ts2, gt2.real, ls='--', color=ax1color, linewidth=linewidth2, label=r'GTEMPO')

ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize, color=ax1color)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize, colors=ax1color)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

# ax[0,0].legend(fontsize=12)

ax2color = 'c'

ax2 = ax[0,0].twinx()



ax2.plot(ts, gt.imag, ls='-', color=ax2color, linewidth=linewidth1, label=r'GTEMPO')
ax2.plot(ts2, gt2.imag, ls='--', color=ax2color, linewidth=linewidth2, label=r'GTEMPO')


ax2.set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))



ax[1,0].plot(ts, lt.real, ls='-', color=ax1color, linewidth=linewidth1, label=r'GTEMPO')
ax[1,0].plot(ts2, lt2.real, ls='--', color=ax1color, linewidth=linewidth2, label=r'GTEMPO')

ax[1,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,0].set_ylabel(r'${\rm Re}[G^{<}(t)]$', fontsize=fontsize, color=ax1color)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize, colors=ax1color)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(c)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax2 = ax[1,0].twinx()

ax2.plot(ts, lt.imag, ls='-', color=ax2color, linewidth=linewidth1, label=r'GTEMPO')
ax2.plot(ts2, lt2.imag, ls='--', color=ax2color, linewidth=linewidth2, label=r'GTEMPO')

ax2.set_xlabel(r'$t$', fontsize=fontsize)
ax2.set_ylabel(r'${\rm Im}[G^{<}(t)]$', fontsize=fontsize, color=ax2color)
ax2.tick_params(axis='both', which='major', labelsize=labelsize, colors=ax2color)
ax2.locator_params(axis='both', nbins=6)
ax2.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

chis = [60,80,100,120,140]

gt_errs = []
lt_errs = []

for chi in chis:
	ts2, taus2, gt2, lt2, gtau2 = read_noninteracting_mixed_tempo(beta, Ntau, t, Nt, mu, d, alpha, chi)
	gt_errs.append(mse_error(gt, gt2))
	lt_errs.append(mse_error(lt, lt2))


ax[0,1].plot(chis, gt_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[0,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$\mathcal{E}[G^{>}(t)]$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)

ax[1,1].plot(chis, lt_errs, ls='--', color='k', linewidth=linewidth2, marker='o', markersize=markersize, markerfacecolor='none')
ax[1,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$\mathcal{E}[G^{<}(t)]$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(d)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

plt.savefig('independentbosons_noint_mixed.pdf', bbox_inches='tight')

plt.show()

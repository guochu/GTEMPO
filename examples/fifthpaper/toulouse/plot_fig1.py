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

def read_mixed_tempo(beta, t, mu, dt, order=10, chi=60):
	filename = 'result/thouless_tempo_mixed_beta%s_dtau0.1_t%s_mu%s_dt%s_order%s_chi%s.json'%(beta, t, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gtau = parse_complex_array(data['gtau'])
	# gf =  gt + lt
	gf = 1j * gt + 1j * lt
	return data['taus'], -gtau.real, data['ts'], gt, lt, gf.imag, data['bd']

def read_ed(beta, N, mu, dt, dw=0.001):
	mpath = '../../../../iGTEMPO/examples/secondpaper/toulouse/'
	filename = mpath + 'result/Toulouse_ed_beta%s_mu%s_dt%s_N%s_dw%s.json'%(beta, mu, dt, N, dw)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	return data['ts'], gt, lt

def read_real_tempo(beta, t, mu, dt, order=10, chi=1024):
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/thouless_tempo_real_beta%s_t%s_mu%s_dt%s_order%s_chi%s.json'%(beta, t, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# gt = parse_complex_array(data['gt'])
	# lt = parse_complex_array(data['lt'])
	# gf =  gt + lt
	return data['ts'], asarray(data['ns']), gt, data['bd']

def read_real_analytic(t, mu, dt):
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/thouless_analytic_real_t%s_mu%s_dt%s.json'%(t, mu, dt)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], asarray(data['gf'])


def read_imag_tempo(beta, mu, dt, order=10, chi=1024):
	N = round(beta/dt)
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/toulouse_tempo_beta%s_mu%s_N%s_order%s_chi%s.json'%(beta, mu, N, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	# lt = asarray(data['lt'])
	return data['ts'], -asarray(data['gf']), asarray(data['bd']).max()

def read_imag_analytic(beta, mu, dt):
	N = round(beta/dt)
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/thouless_analytic_beta%s_mu%s_N%s.json'%(beta, mu, N)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	# lt = asarray(data['lt'])
	return data['ts'], -asarray(data['gf'])

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 14
labelsize = 12
linewidth = 1.5
markersize = 4

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2, 4, figsize=(12,5.5))

t = 40.
beta = 40.
mu = 0.
dt = 0.05
N = round(t / dt)

ts_ed, gt_ed, lt_ed = read_ed(beta, N, mu, dt, dw=0.0005)
ax[0,0].plot(ts_ed, gt_ed.real, ls='--', color='k', linewidth=linewidth, label=r'ED')
ax[0,1].plot(ts_ed, gt_ed.imag, ls='--', color='k', linewidth=linewidth, label=r'ED')


chi = 80
taus, gtau, ts_tempo, gt_tempo, lt_tempo, gf_tempo, bds = read_mixed_tempo(beta, t, mu, dt, chi=chi)
# print(bds)
ax[0,0].plot(ts_tempo, gt_tempo.real, ls='-', color='r', markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi))
ax[0,1].plot(ts_tempo, gt_tempo.imag, ls='-', color='r', markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi))


ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)

ax[0,0].legend(fontsize=10)


ax[0,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,1].set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


chis = [20, 40, 60, 80, 100, 120]

errs_real = []
errs_imag = []
bds = []
diffs = []

for i, chi in enumerate(chis):
	taus, gtau, ts_tempo, gt_tempo, lt_tempo, gf_tempo, bds_mixed = read_mixed_tempo(beta, t, mu, dt, chi=chi)

	err = mse_error(gt_ed.real, gt_tempo.real)
	errs_real.append(err)
	err = mse_error(gt_ed.imag, gt_tempo.imag)
	errs_imag.append(err)

	bds.append(asarray(bds_mixed).max())


ax[1,0].plot(chis, errs_real, ls='--', color='r', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'real')

ax[1,1].plot(chis, errs_imag, ls='--', color='r', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'real')



ax[1,0].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(a2)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)



ax[1,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(b2)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)



chi = 120
taus, gtau, ts_tempo, gt_tempo, lt_tempo, gf_tempo, bds = read_mixed_tempo(beta, t, mu, dt, chi=chi)

ts_real_analytic, gf_real_analytic = read_real_analytic(t, mu, dt)

ax[0,2].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')

ax[0,2].plot(ts_tempo, gf_tempo, ls='--', color=colors[0], linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi))

ax[0,2].set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax[0,2].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,2].locator_params(axis='both', nbins=6)
ax[0,2].annotate(r'(c1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,2].legend(loc='center', fontsize=12)


chis = [10, 20,30, 40,50, 60,70,80]

errs = []
bds = []
diffs = []

for i, chi in enumerate(chis):
	ts_real_tempo, ns_real_tempo, gf_real_tempo, bds_real_tempo = read_real_tempo(beta, t, mu, dt, chi=chi)

	err = mse_error(gf_real_tempo, gf_real_analytic[1:])
	errs.append(err)
	bds.append(asarray(bds_real_tempo).max())

# print(gf_real_analytic)

# print(gf_real_tempo)

ax[1,2].plot(chis, errs, ls='--', color='g', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'real')

chis = [20, 40, 60, 80]
errs = []
bds = []

for i, chi in enumerate(chis):
	taus_mixed, gtau_mixed, ts_mixed, gt_mixed, lt_mixed, gf_mixed, bds_mixed = read_mixed_tempo(beta, t, mu, dt, chi=chi)

	err = mse_error(gf_mixed, gf_real_analytic)
	errs.append(err)
	bds.append(asarray(bds_mixed).max())

ax[1,2].plot(chis, errs, ls='--', color='r', marker=markers[1], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'mixed')


ax[1,2].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,2].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,2].locator_params(axis='both', nbins=6)
ax[1,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,2].annotate(r'(c2)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,2].legend(loc='center', fontsize=12)


# 
chi = 120

dtau = 0.1
ts_imag_analytic, gf_imag_analytic = read_imag_analytic(beta, mu, dt=dtau)

# print(ts_imag_analytic)

ax[0,3].plot(ts_imag_analytic, gf_imag_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')


taus_mixed, gtau_mixed, ts_mixed, gt_mixed, lt_mixed, gf_mixed, bds_mixed = read_mixed_tempo(beta, t, mu, dt=0.05, chi=chi)

ax[0,3].plot(taus_mixed, gtau_mixed, ls='--', color=colors[0], linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi))



ax[0,3].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,3].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,3].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,3].locator_params(axis='both', nbins=6)
ax[0,3].annotate(r'(d1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,3].legend(loc='center', fontsize=12)


chis = [20,40,60,80,100,120]

errs = []
bds = []

for i, chi in enumerate(chis):
	ts_imag_tempo, gf_imag_tempo, bds_real_tempo = read_imag_tempo(beta, mu, dt=dtau, chi=chi)

	err = mse_error(gf_imag_tempo, gf_imag_analytic)
	errs.append(err)
	bds.append(asarray(bds_real_tempo).max())


ax[1,3].plot(bds, errs, ls='--', color='g', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'imag')


chis = [20, 40, 60, 80, 100,120]
errs = []
bds = []

for i, chi in enumerate(chis):
	taus_mixed, gtau_mixed, ts_mixed, gt_mixed, lt_mixed, gf_mixed, bds_mixed = read_mixed_tempo(beta, t, mu, dt, chi=chi)

	err = mse_error(gtau_mixed, gf_imag_analytic)
	errs.append(err)
	bds.append(asarray(bds_mixed).max())

ax[1,3].plot(bds, errs, ls='--', color='r', marker=markers[1], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'mixed')


ax[1,3].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,3].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,3].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,3].locator_params(axis='both', nbins=6)
ax[1,3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax[1,3].annotate(r'(d2)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,3].legend(loc='center', fontsize=12)


plt.tight_layout(pad=0.5)

plt.savefig('toulouse1.pdf', bbox_inches='tight')

plt.show()

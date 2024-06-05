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

def read_mixed_tempo(beta, t, U, dt=0.05, dtau=0.1, order=10, chi=60):
	# mpath = '/Users/guochu/Documents/Missile/GTEMPO/examples/thirdpaper/oneorbital'
	mu = U/2
	# t2 = t/2
	filename = 'result/mixed_beta%s_dtau%s_t%s_dt%s_U%s_mu%s_order%s_chi%s.json'%(beta, dtau, t, dt, U, mu, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gtau = parse_complex_array(data['gtau'])
	gf = 1j * gt + 1j * lt
	ts = asarray(data['ts'])
	return ts-ts[0], gf, gt, lt, data['taus'], -gtau.real


def read_real_tempo(beta, t0, U, dt, order=10, chi=60):
	mpath = '../../thirdpaper/oneorbital/'
	mu = U/2
	# t2 = t/2
	t = t0 + 20.
	filename = mpath + 'result/anderson_tempo1_beta%s_t%s_%s_U%s_mu%s_dt%s_order%s_chi%s_2.json'%(beta, t, t0, U, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = -parse_complex_array(data['lt'])
	gf = 1j * gt + 1j * lt

	gf_ts = dt * asarray([i for i in range(len(data['gf_ts']))])

	return data['ts'], asarray(data['ns']), gf_ts, gf, gt, lt

def read_itempo(beta, t, U, dt, order=6, chi=60, prony=5, k=5):
	N = round(t / dt)
	mu = U/2
	mpath = '/Users/guochu/Documents/Missile/iGTEMPO/examples/secondpaper/oneorbital/'
	filename =  mpath +'result/SIAM_onebath_beta%s_U%s_mu%s_dt%s_k%s_trunc%s_prony%s_N%s_chi%s_it.json'%(beta, U, mu, dt, k, order, prony, N, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gf = 1j * (gt + lt)
	return 0.1 * asarray(data['ts']), gf, gt, lt

def read_imag_tempo(beta, U, dt=0.1, order=10, chi=500):
	mpath = '../../thirdpaper/oneorbital/'
	N = round(beta/dt)
	mu = U/2
	filename = mpath + 'result/anderson_tempo1_beta%s_U%s_mu%s_N%s_dt%s_imag_order%s_chi%s.json'%(beta, U, mu, N, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
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


colors = ['b', 'c', 'orange', 'g', 'y', 'r', 'k']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'r'

fig, ax = plt.subplots(2, 4, figsize=(12,5.5))


beta = 40.
dt = 0.05

U = 0.5

t_final = 80.

t0 = 20.

chi_r = 100

dtau = 0.1

chi_ms = [40,80, 120,160,200]

order = 10

# times_all, ns_final, times_final, gf_final, gt_final, lt_final = read_real_tempo(beta, t_final, U, dt, chi=chi_r)

# ax[0,0].plot(gf_ts_final, gt_final.real, ls='-', color='k', linewidth=linewidth, label=r'real, $\chi=%s$'%(chi_r))
# ax[0,1].plot(gf_ts_final, gt_final.imag, ls='-', color='k', linewidth=linewidth, label=r'real, $\chi=%s$'%(chi_r))

k = 8
times_final, gf_final, gt_final, lt_final = read_itempo(beta, t0, U, dt, chi=chi_r, k=k, order=order)

# ax[0,0].plot(gf_ts_final, gt_final.real, ls='-', color='k', linewidth=linewidth, label=r'real, $\chi=%s$'%(chi_r))
# ax[0,1].plot(gf_ts_final, gt_final.imag, ls='-', color='k', linewidth=linewidth, label=r'real, $\chi=%s$'%(chi_r))


ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


ax[0,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,1].set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)



errs_real = []
errs_imag = []

errs = []

for i, chi_m in enumerate(chi_ms):
	

	mixed_ts, mixed_gf, mixed_gt, mixed_lt, mixed_taus, mixed_gtau = read_mixed_tempo(beta, t0, U, dt, dtau=dtau, chi=chi_m)

	ax[0,0].plot(mixed_ts, mixed_gt.real, ls='--', color=colors[i], linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi_m))

	ax[0,1].plot(mixed_ts, mixed_gt.imag, ls='--', color=colors[i], linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi_m))

	err = mse_error(mixed_gt.real, gt_final.real)

	errs_real.append(err)

	err = mse_error(mixed_gt.imag, gt_final.imag)

	errs_imag.append(err)

	err = mse_error(mixed_gf, gf_final)

	errs.append(err)

ax[0,0].legend(fontsize=10)

ax[1,0].plot(chi_ms, errs_real, ls='--', color='r', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'real')

ax[1,1].plot(chi_ms, errs_imag, ls='--', color='r', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'real')


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


chi = 100
dtau = 0.1

for i, chi_m in enumerate(chi_ms):
	mixed_ts, mixed_gf, mixed_gt, mixed_lt, mixed_taus, mixed_gtau = read_mixed_tempo(beta, t0, U, dt, dtau=dtau, chi=chi_m)
	ax[0,2].plot(mixed_ts, gf_final.imag, ls='-', color='k', linewidth=linewidth)
	ax[0,2].plot(mixed_ts, mixed_gf.imag, ls='--', color=colors[i], linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi_m))


ax[0,2].set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax[0,2].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,2].locator_params(axis='both', nbins=6)
ax[0,2].annotate(r'(c1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,2].legend(fontsize=10)


ax[1,2].plot(chi_ms, errs, ls='--', color='r', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'real')

ax[1,2].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,2].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,2].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,2].locator_params(axis='both', nbins=6)
ax[1,2].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,2].annotate(r'(c2)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


chi_i = 500

taus, gtau = read_imag_tempo(beta, U)
ax[0, 3].plot(taus, gtau, ls='-', color='k', linewidth=linewidth, label=r'imag, $\chi=%s$'%(chi_i))


errs_i = []

for i, chi_m in enumerate(chi_ms):
	

	mixed_ts, mixed_gf, mixed_gt, mixed_lt, mixed_taus, mixed_gtau = read_mixed_tempo(beta, t0, U, dt, chi=chi_m)

	ax[0,3].plot(mixed_taus, mixed_gtau, ls='--', color=colors[i], linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi_m))


	err = mse_error(mixed_gtau, gtau)

	errs_i.append(err)



ax[0,3].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,3].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,3].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,3].locator_params(axis='both', nbins=6)
ax[0,3].annotate(r'(d1)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,3].legend(loc='center', fontsize=10)


ax[1,3].plot(chi_ms, errs_i, ls='--', color='r', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none', label=r'mixed')


ax[1,3].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,3].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,3].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,3].locator_params(axis='both', nbins=6)
ax[1,3].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax[1,3].annotate(r'(d2)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[1,3].legend(loc='center', fontsize=12)


plt.tight_layout(pad=0.5)

# plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


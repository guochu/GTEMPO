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

def read_mixed_tempo(beta, t, U, dt=0.05, order=10, chi=60):
	# mpath = '/Users/guochu/Documents/Missile/GTEMPO/examples/thirdpaper/oneorbital'
	mu = U/2
	# t2 = t/2
	filename = 'result/mixed_beta%s_dtau0.1_t%s_dt%s_U%s_mu%s_order%s_chi%s.json'%(beta, t, dt, U, mu, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gf = 1j * gt + 1j * lt
	ts = asarray(data['ts'])
	return ts-ts[0], gf, gt, lt


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


def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)


fontsize = 20
labelsize = 16
linewidth = 2.
markersize = 10

colors = ['b', 'c', 'orange', 'g', 'y', 'r', 'k']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'r'

fig, ax = plt.subplots(1, 2, figsize=(8, 4))


beta = 40.
dt = 0.05

U = 0.5

t_final = 80.

t0 = 20.

chi_r = 60


chi_ms = [40, 60, 80, 100, 120, 140]

times_final, ns_final, gf_ts_final, gf_final, gt_final, lt_final = read_real_tempo(beta, t_final, U, dt, chi=chi_r)

ax[0].plot(gf_ts_final, gf_final.imag, ls='-', color='k', linewidth=linewidth, label=r'real, $\chi=%s$'%(chi_r))


errs = []

for i, chi_m in enumerate(chi_ms):
	

	mixed_ts, mixed_gf, mixed_gt, mixed_lt = read_mixed_tempo(beta, t0, U, dt, chi=chi_m)

	ax[0].plot(mixed_ts, mixed_gf.imag, ls='--', color=colors[i], linewidth=linewidth, label=r'imag, $\chi=%s$'%(chi_m))

	errs.append(mse_error(gf_final.imag, mixed_gf.imag))

	print(mse_error(gf_final.imag, mixed_gf.imag))


ax[0].set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax[0].set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].tick_params(axis='y', which='both')
ax[0].locator_params(axis='both', nbins=6)


ax[0].legend(fontsize=12)



ax[1].plot(chi_ms, errs, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth)

ax[1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].tick_params(axis='y', which='both')
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1].locator_params(axis='both', nbins=6)

plt.tight_layout(pad=0.5)

# plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


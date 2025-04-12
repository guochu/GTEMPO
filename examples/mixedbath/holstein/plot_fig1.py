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


def read_real_tempo(t, N, mu, omega, alpha, chi=80):
	dt = t / N
	filename = 'result/holstein_realgtempo_betaInf_t%s_dt%s_omega0%s_alpha0%s_mu%s_chi%s.json'%(t, dt, omega, alpha, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	ts = asarray(data['ts'])
	return ts-ts[0], -gt * 1j, data['bd1'], data['bd2']


def read_real_analytic(t, N, mu, omega, alpha, order=10, wmax=100.):
	dt = t / N
	filename = 'result/holstein_analytic_betaInf_t%s_dt%s_omega0%s_alpha0%s_mu%s_order%s_wmax%s.json'%(t, dt, omega, alpha, mu, order, wmax)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gf = parse_complex_array(data['gf'])
	ts = asarray(data['ts'])
	return ts-ts[0], gf 



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

fig, ax = plt.subplots(1,2, figsize=(8,3.5))


chi = 100

mu = -1.1

t = 5
Nt = 400
omega = 1
alpha = 0.5

order = 20
wmax = 50

# noninteracting case
ts, gt, bd1, bd2 = read_real_tempo(t, Nt, mu, omega, alpha, chi)

ax[0].plot(ts, gt.real, ls='--', color='k', linewidth=linewidth, label=r'GTEMPO')
ax[1].plot(ts, gt.imag, ls='--', color='k', linewidth=linewidth, label=r'GTEMPO')


Nt = 400
ts2, gt2 = read_real_analytic(t, Nt, mu, omega, alpha, order=order, wmax=wmax)

ax[0].plot(ts2, gt2.real, ls='-', color='r', linewidth=linewidth, label=r'Analytic')
ax[1].plot(ts2, gt2.imag, ls='-', color='r', linewidth=linewidth, label=r'Analytic')

# print(len(gt), len(gt2))
# print('errors: ', mse_error(gt, gt2))


ax[0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0].set_ylabel(r'${\rm Re}[G(t)]$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].legend(loc='lower right', fontsize=12)

ax[0].legend(fontsize=12)


ax[1].set_xlabel(r'$t$', fontsize=fontsize)
ax[1].set_ylabel(r'${\rm Im}[G(t)]$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].locator_params(axis='both', nbins=6)
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

# plt.savefig('fig2.pdf', bbox_inches='tight')

plt.show()

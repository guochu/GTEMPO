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

def read_interacting_imag_tempo(beta, Ntau, U, mu, d, alpha, chi=80):
	dtau = beta / Ntau
	filename = 'result/interacting_imaggtempo_beta%s_dtau%s_U%s_mu%s_d%s_alpha%s_chi%s.json'%(beta, dtau, U, mu, d, alpha, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gtau = asarray(data['gtau'])
	taus = asarray(data['taus'])
	return taus-taus[0], gtau

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

fontsize = 20
labelsize = 18
linewidth = 2.5
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(1,1, figsize=(8,6))


chi = 160

U = 1
mu = 0.5

beta = 5
Ntau = 25
d = 1
alpha = 1

# noninteracting case
taus, gtau = read_interacting_imag_analytic(beta, Ntau, U, mu, d, alpha)

taus2, gtau2 = read_interacting_imag_tempo(beta, Ntau, U, mu, d, alpha, chi)
# ts4, taus4, gt4, lt4, gtau4 = read_noninteracting_mixed_tempo(beta, Ntau, 0.1, 10, mu, d, alpha, chi)

ax.plot(taus, gtau, ls='-', color='k', linewidth=linewidth, label=r'Analytic')
ax.plot(taus2, gtau2.real, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')

ax.set_xlabel(r'$\tau$', fontsize=fontsize)
ax.set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.locator_params(axis='both', nbins=6)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

gtau_errors = []
chis = [20, 40,60,80,100, 120, 140, 160]

for i, chi in enumerate(chis):
	taus2, gtau2 = read_interacting_imag_tempo(beta, Ntau, U, mu, d, alpha, chi)
	gtau_errors.append(mse_error(gtau, gtau2.real))


linewidth_s = 1.4
markersize_s = 4
fontsize_s = 16
labelsize_s = 14

ax1 = ax.inset_axes([0.25, 0.2, 0.5, 0.5])

ax1.plot(chis, gtau_errors, ls='--', color='k', marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth_s, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\chi$', fontsize=fontsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))



ax.legend(fontsize=14)

plt.tight_layout(pad=0.5)

plt.show()

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

fontsize = 14
labelsize = 12
linewidth = 1.5
markersize = 4

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2,2, figsize=(8,6))


chi = 100

mu = 0.

t = 2
Nt = 200
beta = 0.5
Ntau = 50
d = 1
alpha = 1

# noninteracting case
taus, gtau = read_noninteracting_imag_analytic(beta, Ntau, mu, d, alpha)
ts, gt, lt = read_noninteracting_real_analytic(beta, t, Nt, mu, d, alpha)
# gf = gt - lt

ax[0,0].plot(taus, gtau, ls='--', color='k', linewidth=linewidth, label=r'Analytic')
ax[0,1].plot(ts, lt.imag, ls='--', color='k', linewidth=linewidth, label=r'Analytic')


ts2, taus2, gt2, lt2, gtau2 = read_noninteracting_mixed_tempo(beta, Ntau, t, Nt, mu, d, alpha, chi)
# gf2 = gt2 - lt2

taus3, gtau3 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)

ts4, taus4, gt4, lt4, gtau4 = read_noninteracting_mixed_tempo(beta, Ntau, 0.1, 10, mu, d, alpha, chi)


ax[0,0].plot(taus2, gtau2.real, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')
ax[0,1].plot(ts, lt2.imag, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')

print('errors: ', mse_error(gtau, gtau2), ' ', mse_error(gt, gt2))

print(gtau[:10])
print(gtau2[:10])
# print(asarray(gtau) - gtau2 / (gtau2[0] + gtau2[-1]))

ax[0,0].legend(fontsize=12)


# # interacting case
# U = 1

# taus, gtau = read_interacting_imag_analytic(beta, Ntau, U, mu, d, alpha)
# ts, gt, lt = read_interacting_real_analytic(beta, t, Nt, U, mu, d, alpha)
# gf = gt - lt

# ax[1,0].plot(taus, gtau, ls='--', color='k', linewidth=linewidth, label=r'Analytic')
# ax[1,1].plot(ts, lt.imag, ls='--', color='k', linewidth=linewidth, label=r'Analytic')


# ts2, taus2, gt2, lt2, gtau2 = read_interacting_mixed_tempo(beta, Ntau, t, Nt, U, mu, d, alpha, chi)
# gf2 = gt2 - lt2

# ax[1,0].plot(taus2, gtau2.real, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')
# ax[1,1].plot(ts, lt2.imag, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')

# print('errors: ', mse_error(gtau, gtau2), ' ', mse_error(gt, gt2))

ax[1,0].legend(fontsize=12)


plt.show()

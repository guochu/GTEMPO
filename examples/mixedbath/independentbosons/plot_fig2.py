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


def read_noninteracting_mixed_tempo(beta, Ntau, t, N, mu, omega, alpha, chi=80):
	dt = t / N
	dtau = beta / Ntau
	filename = 'result/noninteracting_sm_mixedgtempo_beta%s_dtau%s_t%s_dt%s_mu%s_omega%s_alpha%s_chi%s.json'%(beta, dtau, t, dt, mu, omega, alpha, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gtau = parse_complex_array(data['gtau'])
	ts = asarray(data['ts'])
	taus = asarray(data['taus'])
	return ts-ts[0], taus-taus[0], gt, lt, gtau

def read_noninteracting_imag_analytic(beta, N, mu, omega, alpha):
	filename = 'result/noninteracting_sm_analytic_imag_beta%s_mu%s_N%s_omega%s_alpha%s.json'%(beta, mu, N, omega, alpha)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = data['gf']
	return data['ts'], gt

def read_noninteracting_real_analytic(beta, t, N, mu, omega, alpha):
	filename = 'result/noninteracting_sm_analytic_real_beta%s_mu%s_t%s_N%s_omega%s_alpha%s.json'%(beta, mu, t, N, omega, alpha)
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
omega = 1
alpha = 0.5

# noninteracting case
taus, gtau = read_noninteracting_imag_analytic(beta, Ntau, mu, omega, alpha)
ts, gt, lt = read_noninteracting_real_analytic(beta, t, Nt, mu, omega, alpha)
# gf = gt - lt

ax[0,0].plot(taus, gtau, ls='--', color='k', linewidth=linewidth, label=r'Analytic')
ax[0,1].plot(ts, lt.imag, ls='--', color='k', linewidth=linewidth, label=r'Analytic')


ts2, taus2, gt2, lt2, gtau2 = read_noninteracting_mixed_tempo(beta, Ntau, t, Nt, mu, omega, alpha, chi)
# gf2 = gt2 - lt2

ax[0,0].plot(taus2, gtau2.real, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')
ax[0,1].plot(ts, lt2.imag, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')

print('errors: ', mse_error(gtau, gtau2), ' ', mse_error(gt, gt2))

ax[0,0].legend(fontsize=12)


ax[1,0].legend(fontsize=12)


plt.show()

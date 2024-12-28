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
	return ts-ts[0], gt

# def read_noninteracting_imag_analytic(beta, N, mu, d, alpha):
# 	filename = 'result/noninteracting_analytic_imag_beta%s_mu%s_N%s_d%s_alpha%s.json'%(beta, mu, N, d, alpha)
# 	with open(filename, 'r') as f:
# 		data = f.read()
# 		data = json.loads(data)
# 	gt = data['gf']
# 	return data['ts'], gt

# def read_noninteracting_real_analytic(beta, t, N, mu, d, alpha):
# 	filename = 'result/noninteracting_analytic_real_beta%s_mu%s_t%s_N%s_d%s_alpha%s.json'%(beta, mu, t, N, d, alpha)
# 	with open(filename, 'r') as f:
# 		data = f.read()
# 		data = json.loads(data)
# 	gt = parse_complex_array(data['gt'])
# 	lt = parse_complex_array(data['lt'])
# 	return data['ts'], gt, lt



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
omega = 1
alpha = 1

# noninteracting case
ts, gt = read_real_tempo(t, Nt, mu, omega, alpha, chi)

ax[0,0].plot(ts, gt.real, ls='--', color='k', linewidth=linewidth, label=r'GTEMPO')

ax[0,1].plot(ts, gt.imag, ls='--', color='k', linewidth=linewidth, label=r'GTEMPO')


ax[0,0].legend(fontsize=12)



ax[1,0].legend(fontsize=12)


plt.show()

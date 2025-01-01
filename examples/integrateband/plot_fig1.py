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

def read_tempo_1(beta, N, U, mu, chi=80):
	filename = "result/anderson_tempo_beta%s_U%s_e%s_N%s_chi%s.json"%(beta, U, mu, N, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gtau = asarray(data['gf'])
	ts = asarray(data['ts'])
	return ts, gtau


def read_tempo_2(beta, N, U, mu, chi=80, chi2=500):
	filename = "result/anderson_tempo_beta%s_U%s_e%s_N%s_chi%s_chi2%s.json"%(beta, U, mu, N, chi, chi2)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gtau = asarray(data['gf'])
	ts = asarray(data['ts'])
	return ts, gtau


def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 20
labelsize = 18
linewidth = 2.5
markersize = 4

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(1,1, figsize=(8,6))



beta = 20.
dtau = 0.1
N = round(beta/dtau)

U = 1.
mu = U/2
chi = 60
chi2 = 500
# Ntau = 20

# real time data
ts, gtau = read_tempo_1(beta, N, U, mu, chi)
ax.plot(ts, gtau, ls='-', color='k', linewidth=linewidth, label=r'GTEMPO1')

ts2, gtau2 = read_tempo_2(beta, N, U, mu, chi, chi2)
ax.plot(ts2, gtau2, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO1')

print(mse_error(gtau, gtau2))

ax.set_xlabel(r'$\tau$', fontsize=fontsize)
ax.set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.locator_params(axis='both', nbins=6)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax.legend(loc='lower right', fontsize=12)




plt.tight_layout(pad=0.5)

plt.show()

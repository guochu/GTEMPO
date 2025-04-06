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

def read_real_tempo(beta, t, N, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dt = t / N
	filename = 'result/noninteracting_realgtempo_beta%s_t%s_dt%s_omega0%s_alpha0%s_omega1%s_alpha1%s_mu%s_chi%s.json'%(beta, t, dt, omega0, alpha0, omega1, alpha1, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt

def read_mixed_tempo(beta, Ntau, t, N, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dt = t / N
	dtau = beta / Ntau
	filename = 'result/noninteracting_mixedgtempo_beta%s_dtau%s_t%s_dt%s_omega0%s_alpha0%s_omega1%s_alpha1%s_mu%s_chi%s.json'%(beta, dtau, t, dt, omega0, alpha0, omega1, alpha1, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt

def read_neq_ed(beta, t, N, mu, omega0, alpha0, omega1, alpha1):
	filename = 'result/noninteracting_neq_ED_real_beta%s_mu%s_t%s_N%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, t, N, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	return data['ts'], gt, lt

def read_eq_ed(beta, t, N, mu, omega0, alpha0, omega1, alpha1):
	filename = 'result/noninteracting_eq_ED_real_beta%s_mu%s_t%s_N%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, t, N, omega0, alpha0, omega1, alpha1)
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


omega0 = 1
alpha0 = 0.5
omega1 = 1
alpha1 = 1
chi = 100

mu = 0.

t = 10
Nt = 200
beta = 5
# Ntau = 20

# real time data
ts, gt, lt = read_neq_ed(beta, t, Nt, mu, omega0, alpha0, omega1, alpha1)

ax[0,0].plot(ts, gt.real, ls='--', color='k', linewidth=linewidth, label=r'ED')
ax[0,1].plot(ts, gt.imag, ls='--', color='k', linewidth=linewidth, label=r'ED')

ax[1,0].plot(ts, lt.real, ls='--', color='k', linewidth=linewidth, label=r'ED')
ax[1,1].plot(ts, lt.imag, ls='--', color='k', linewidth=linewidth, label=r'ED')

ts2, gt2, lt2 = read_real_tempo(beta, t, Nt, mu, omega0, alpha0, omega1, alpha1, chi)

print('error is ', mse_error(gt, gt2), ' ', mse_error(lt, lt2))

ax[0,0].plot(ts2, gt2.real, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')
ax[0,1].plot(ts2, gt2.imag, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')


ax[1,0].plot(ts2, lt2.real, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')
ax[1,1].plot(ts2, lt2.imag, ls='--', color='r', linewidth=linewidth, label=r'GTEMPO')

ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].set_ylabel(r'${\rm Re}[G^{>}(t)]$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(loc='lower right', fontsize=12)


ax[0,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,1].set_ylabel(r'${\rm Im}[G^{>}(t)]$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)

# 
ax[1,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,0].set_ylabel(r'${\rm Re}[G^{<}(t)]$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(c)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[1,0].legend(fontsize=12)


ax[1,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,1].set_ylabel(r'${\rm Im}[G^{<}(t)]$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(d)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)



plt.tight_layout(pad=0.5)

plt.show()

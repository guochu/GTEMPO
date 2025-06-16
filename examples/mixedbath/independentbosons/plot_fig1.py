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

def read_noninteracting_imag_tempo_2(beta, Ntau, mu, d, alpha, chi=80):
	dtau = beta / Ntau
	filename = 'result/noninteracting_imaggtempo_beta%s_dtau%s_mu%s_d%s_alpha%s_chi%s_2.json'%(beta, dtau, mu, d, alpha, chi)
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

fontsize = 22
labelsize = 18
linewidth1 = 1.5
linewidth2 = 3.
markersize = 10

linewidth_s = 2
markersize_s = 4
fontsize_s = 16
labelsize_s = 14


colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2,2, figsize=(8,7))

annotate_xy = (-0.15, 1.05)

chi = 140


beta = 1
Ntau = 20
d = 1
alpha = 1

mu = 0.

color = 'k'
color2 = 'g'

# noninteracting case
taus, gtau = read_noninteracting_imag_analytic(beta, Ntau, mu, d, alpha)

taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)


ax[0,0].plot(taus, gtau, ls='-', color=color, linewidth=linewidth1, label=r'Analytic')
ax[0,0].plot(taus2, gtau2.real, ls='--', color=color, linewidth=linewidth2, label=r'GTEMPO')


ax[0,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].set_title(r'$\epsilon_d=%s, d=%s, \alpha=%s$'%(round(-mu), d, alpha), fontsize=labelsize_s)
ax[0,0].annotate(r'(a)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


gtau_errors = []
chis = [20, 40,60,80, 100, 120, 140]

for i, chi in enumerate(chis):
	taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)
	gtau_errors.append(mse_error(gtau, gtau2.real))


ax1 = ax[0,0].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.semilogy(chis, gtau_errors, ls='--', color=color, marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')


ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\chi$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)




# 
beta = 10
Ntau = 50

taus, gtau = read_noninteracting_imag_analytic(beta, Ntau, mu, d, alpha)

taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)


ax[0,1].plot(taus, gtau, ls='-', color=color, linewidth=linewidth1, label=r'Analytic')
ax[0,1].plot(taus2, gtau2.real, ls='--', color=color, linewidth=linewidth2, label=r'GTEMPO')

ax[0,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].set_title(r'$\epsilon_d=%s, d=%s, \alpha=%s$'%( round(-mu), d, alpha), fontsize=labelsize_s)

ax[0,1].annotate(r'(b)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


gtau_errors = []
# chis = [20, 40,60,80, 100, 120, 140]

for i, chi in enumerate(chis):
	taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)
	gtau_errors.append(mse_error(gtau, gtau2.real))



ax1 = ax[0,1].inset_axes([0.25, 0.4, 0.5, 0.5])

ax1.semilogy(chis, gtau_errors, ls='--', color=color, marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\chi$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)


# 
mu = -5.

taus, gtau = read_noninteracting_imag_analytic(beta, Ntau, mu, d, alpha)

taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)

ax[1,0].plot(taus, gtau, ls='-', color=color, linewidth=linewidth1, label=r'Analytic')
ax[1,0].plot(taus2, gtau2.real, ls='--', color=color, linewidth=linewidth2, label=r'GTEMPO')

ax[1,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].set_title(r'$\epsilon_d=%s, d=%s, \alpha=%s$'%( round(-mu), d, alpha), fontsize=labelsize_s)

ax[1,0].annotate(r'(c)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


gtau_errors = []
# chis = [20, 40,60,80, 100, 120, 140]

for i, chi in enumerate(chis):
	taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)
	gtau_errors.append(mse_error(gtau, gtau2.real))



ax1 = ax[1,0].inset_axes([0.3, 0.4, 0.5, 0.5])

ax1.semilogy(chis, gtau_errors, ls='--', color=color, marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\chi$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)



# 
d = 3

taus, gtau = read_noninteracting_imag_analytic(beta, Ntau, mu, d, alpha)

taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)

ax[1,1].plot(taus, gtau, ls='-', color=color, linewidth=linewidth1, label=r'Analytic')
ax[1,1].plot(taus2, gtau2.real, ls='--', color=color, linewidth=linewidth2, label=r'GTEMPO')

ax[1,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].set_title(r'$\epsilon_d=%s, d=%s, \alpha=%s$'%( round(-mu), d, alpha), fontsize=labelsize_s)

ax[1,1].annotate(r'(d)', xy=annotate_xy,xycoords='axes fraction', fontsize=fontsize)


gtau_errors = []
# chis = [20, 40,60,80, 100, 120, 140]

for i, chi in enumerate(chis):
	taus2, gtau2 = read_noninteracting_imag_tempo(beta, Ntau, mu, d, alpha, chi)
	gtau_errors.append(mse_error(gtau, gtau2.real))



ax1 = ax[1,1].inset_axes([0.3, 0.4, 0.5, 0.5])

ax1.semilogy(chis, gtau_errors, ls='--', color=color, marker='o', markersize=markersize, markerfacecolor='none', linewidth=linewidth1, label=r'Partial')

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize_s)
ax1.set_xlabel(r'$\chi$', fontsize=fontsize_s)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)



# ax[0,0].legend(fontsize=14)

plt.tight_layout(pad=0.5)

plt.savefig('independentbosons_noint_imag.pdf', bbox_inches='tight')

plt.show()

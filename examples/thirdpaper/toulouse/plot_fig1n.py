import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace
from numpy.linalg import norm
import math


rc('text', usetex=True)


def read_real_tempo(beta, t, mu, dt, order=10, chi=1024):
	filename = 'result/thouless_tempo_real_beta%s_t%s_mu%s_dt%s_order%s_chi%s.json'%(beta, t, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# lt = asarray(data['lt'])
	return data['ts'], asarray(data['ns']), gt, data['bd']

def read_real_analytic(t, mu, dt):
	filename = 'result/thouless_analytic_real_t%s_mu%s_dt%s.json'%(t, mu, dt)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], asarray(data['gf'])


def read_imag_tempo(beta, mu, dt, order=10, chi=1024):
	N = round(beta/dt)
	filename = 'result/toulouse_tempo_beta%s_mu%s_N%s_order%s_chi%s.json'%(beta, mu, N, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	# lt = asarray(data['lt'])
	return data['ts'], -asarray(data['gf']), asarray(data['bd']).max()

def read_imag_analytic(beta, mu, dt):
	N = round(beta/dt)
	filename = 'result/thouless_analytic_beta%s_mu%s_N%s.json'%(beta, mu, N)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	# lt = asarray(data['lt'])
	return data['ts'], -asarray(data['gf'])

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return v * v / L

ts_real_analytic, gf_real_analytic = read_real_analytic(40., 1., 0.05)

ts, ns, gf, bds = read_real_tempo(40., 40., 1., 0.05, 7)

# print(gf_real_analytic[-20:])
# print(gf[-20:])

# ts, gf, bds = read_imag_tempo(40., 0., 0.1, 8)

# print(gf)

fontsize = 20
labelsize = 16
linewidth = 1.5
markersize = 7

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2, 2, figsize=(8, 7))

t = 40.
beta = 40.
mu = 0.
dt = 0.05

ts_real_analytic, gf_real_analytic = read_real_analytic(t, mu, dt)

ax[0,0].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth)

chis = [10, 20, 30, 40]

errs = []
bds = []
diffs = []

for i, chi in enumerate(chis):
	ts_real_tempo, ns_real_tempo, gf_real_tempo, bds_real_tempo = read_real_tempo(beta, t, mu, dt, chi=chi)
	ax[0,0].plot(ts_real_tempo, gf_real_tempo, ls='--', color=colors[i], linewidth=linewidth, label=r'$\chi=%s$'%(chi))

	err = mse_error(gf_real_tempo, gf_real_analytic[1:])
	errs.append(err)
	diffs.append(gf_real_analytic[1:] - gf_real_tempo)
	bds.append(asarray(bds_real_tempo).max())



ax[0,0].set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(loc='center', fontsize=12)

# ax[0,0].set_ylim(top=1.2)
# ax[0,0].legend(loc='upper center', fontsize=12)


ax[0,1].semilogy(chis, errs, ls='--', color='k', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none')


ax[0,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[0,1].set_ylabel(r'MSE error', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[0,1].locator_params(axis='both', nbins=6)
# ax[0,1].legend(loc='center right', fontsize=12)
ax[0,1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)



# ax1 = ax[0,0].inset_axes([0.35, 0.4, 0.45, 0.45])

# linewidth_s = 1.
# markersize_s = 4
# fontsize_s = 16
# labelsize_s = 12


# for (i, chi) in enumerate(chis):
# 	ax1.plot(ts_real_analytic[1:], diffs[i], linewidth=linewidth_s, color=colors[i], markersize=markersize_s, markerfacecolor='none', ls='--', label=r'$\chi=%s$'%(chi))


# ax1.set_ylabel(r'Error', fontsize=fontsize_s)
# ax1.set_xlabel(r'$t$', fontsize=fontsize_s, labelpad=0.1)
# ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax1.locator_params(axis='both', nbins=6)


chis = [20,40,60,80,100]

mu = 0.
dt = 0.1
ts_imag_analytic, gf_imag_analytic = read_imag_analytic(beta, mu, dt)

# print(ts_imag_analytic)

ax[1,0].plot(ts_imag_analytic, gf_imag_analytic, ls='-', color='k', linewidth=linewidth)



for i, chi in enumerate(chis):
	ts_imag_tempo, gf_imag_tempo, bds_real_tempo = read_imag_tempo(beta, mu, dt, chi=chi)

	ax[1,0].plot(ts_imag_tempo, gf_imag_tempo, ls='--', color=colors[i], linewidth=linewidth, label=r'$\chi=%s$'%(chi))



ax[1,0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].annotate(r'(c)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,0].legend(loc='center', fontsize=12)


chis = [20,40,60,80,100, 120, 140, 160]

errs = []
bds = []

for i, chi in enumerate(chis):
	ts_imag_tempo, gf_imag_tempo, bds_real_tempo = read_imag_tempo(beta, mu, dt, chi=chi)

	err = mse_error(gf_imag_tempo, gf_imag_analytic)
	errs.append(err)
	bds.append(asarray(bds_real_tempo).max())


ax[1,1].semilogy(bds, errs, ls='--', color='k', marker=markers[0], markersize=markersize, linewidth=linewidth, markerfacecolor='none')


ax[1,1].set_xlabel(r'$\chi$', fontsize=fontsize)
ax[1,1].set_ylabel(r'MSE error', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[0,1].locator_params(axis='both', nbins=6)
# ax[0,1].legend(loc='center right', fontsize=12)
ax[1,1].annotate(r'(d)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)




# ax1 = ax[1,0].inset_axes([0.3, 0.2, 0.45, 0.45])

# linewidth_s = 1.
# markersize_s = 4
# fontsize_s = 16
# labelsize_s = 12

# # diff1 = abs(gf - gf_tempo)
# diff_e8 = gf_imag_analytic - gf_imag_tempo_e8
# diff_e9 = gf_imag_analytic - gf_imag_tempo_e9
# diff_e10 = gf_imag_analytic - gf_imag_tempo_e10

# # ax1.plot(ts_imag_tempo_e8, diff_e8, linewidth=linewidth_s, color=colors[2], markersize=markersize_s, markerfacecolor='none', ls='--')
# ax1.plot(ts_imag_tempo_e9, diff_e9, linewidth=linewidth_s, color=colors[3], markersize=markersize_s, markerfacecolor='none', ls='--', label=r'$\varsigma=10^{-9}$')
# ax1.plot(ts_imag_tempo_e10, diff_e10, linewidth=linewidth_s, color=colors[4], markersize=markersize_s, markerfacecolor='none', ls='--', label=r'$\varsigma=10^{-10}$')


# ax1.set_ylabel(r'Error', fontsize=fontsize_s)
# ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
# ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax1.locator_params(axis='both', nbins=6)


# betas = [5., 10., 15., 20., 25., 30., 35., 40.]
# bds_e8 = [read_imag_tempo(beta, mu, dt, 8)[2] for beta in betas]
# bds_e9 = [read_imag_tempo(beta, mu, dt, 9)[2] for beta in betas]
# bds_e10 = [read_imag_tempo(beta, mu, dt, 10)[2] for beta in betas]

# ax[1,1].plot(betas, bds_e8, ls='--', color=colors[2], linewidth=linewidth, marker=markers[0], markersize=markersize, markerfacecolor='none', label=r'$\varsigma=10^{-8}$')
# ax[1,1].plot(betas, bds_e9, ls='--', color=colors[3], linewidth=linewidth, marker=markers[1], markersize=markersize, markerfacecolor='none', label=r'$\varsigma=10^{-9}$')
# ax[1,1].plot(betas, bds_e10, ls='--', color=colors[4], linewidth=linewidth, marker=markers[2], markersize=markersize, markerfacecolor='none', label=r'$\varsigma=10^{-10}$')


# ax[1,1].set_ylabel(r'$\chi$', fontsize=fontsize)
# ax[1,1].set_xlabel(r'$\beta$', fontsize=fontsize)
# ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[1,1].locator_params(axis='both', nbins=6)

# # ax[1,1].legend(loc='center right', fontsize=12)
# ax[1,1].annotate(r'(d)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)



plt.tight_layout(pad=0.5)

# plt.savefig('toulouse1.pdf', bbox_inches='tight')

plt.show()





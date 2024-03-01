import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace
import math


rc('text', usetex=True)

def parse_complex_array(data):
	re = [item['re'] for item in data]
	im = [item['im'] for item in data]
	return asarray(re) + 1j * asarray(im)


def read_real_tempo(beta, t, mu, order=10, chi=60):
	t2 = 100.
	lb = -2.
	ub = 2.
	dt = 0.05
	filename = 'result/thouless_tempo_real_beta%s_t%s_mu%s_tep%s_lb%s_ub%s_dt0.05_order%s_chi%s.json'%(beta, t, mu, t2, lb, ub, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)

	return data['ts'], data['gf'], data['ws'], data['Aw']

def read_imag_tempo(beta, mu, dt=0.1, order=10, chi=1024):
	N = round(beta/dt)
	filename = 'result/toulouse_tempo_beta%s_mu%s_N%s_order%s_chi%s.json'%(beta, mu, N, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], -asarray(data['gf'])

def read_imag_analytic(beta, mu, dt=0.1):
	N = round(beta/dt)
	filename = 'result/thouless_analytic_beta%s_mu%s_N%s.json'%(beta, mu, N)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	# lt = asarray(data['lt'])
	return data['ts'], -asarray(data['gf'])





fontsize = 20
labelsize = 16
linewidth = 1.5
markersize = 7

colors = ['b', 'c', 'y']
markers = ['o', '^', '+']

color0 = 'r'
color1 = 'g'
color2 = 'tab:orange'
marker1 = 'o'
marker2 = '+'


fig, ax = plt.subplots(1, 2, figsize=(8, 3.5))

beta = 40.
t = 40.
mu = 0.



ts_imag_analytic, gf_imag_analytic = read_imag_analytic(beta, mu)



ax[0].plot(ts_imag_analytic, gf_imag_analytic, ls='-', color='k', linewidth=linewidth)

chi = 40
ts, gf, ws, Aw = read_real_tempo(beta, t, mu, chi=chi)

ax[0].plot(ts, gf, ls='--', color=color1, linewidth=linewidth, label=r'real, $\chi=%s$'%(chi))

chi = 120
ts_imag, gf_imag = read_imag_tempo(beta, mu, chi=chi)

ax[0].plot(ts_imag, gf_imag, ls='--', color=color0, linewidth=linewidth, label=r'imag, $\chi=%s$'%(chi))


ax[0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
ax[0].annotate(r'(a)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].legend(loc='upper right', fontsize=12)

ax1 = ax[0].inset_axes([0.35, 0.2, 0.45, 0.45])

linewidth_s = 1.
markersize_s = 4
fontsize_s = 16
labelsize_s = 12


diff = gf_imag_analytic - gf

diff_imag = gf_imag_analytic - gf_imag

ax1.plot(ts, diff, linewidth=linewidth_s, color=color1, markersize=markersize_s, markerfacecolor='none', ls='--')

ax1.plot(ts_imag, diff_imag, linewidth=linewidth_s, color=color0, markersize=markersize_s, markerfacecolor='none', ls='--')


ax1.set_ylabel(r'Error', fontsize=fontsize_s)
ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax1.locator_params(axis='both', nbins=5)



mu = 1.

# ts_imag_analytic, gf_imag_analytic = read_imag_analytic(beta, mu)



# ax[1].plot(ts_imag_analytic, gf_imag_analytic, ls='-', color='k', linewidth=linewidth)

# chi = 40
# ts, gf, ws, Aw = read_real_tempo(beta, t, mu, chi=chi)

# ax[1].plot(ts, gf, ls='--', color=color1, linewidth=linewidth, label=r'real, $\chi=%s$'%(chi))

# chi = 120
# ts_imag, gf_imag = read_imag_tempo(beta, mu, chi=chi)

# ax[1].plot(ts_imag, gf_imag, ls='--', color=color0, linewidth=linewidth, label=r'imag, $\chi=%s$'%(chi))


# ax[1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
# ax[1].set_xlabel(r'$\tau$', fontsize=fontsize)
# ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[1].locator_params(axis='both', nbins=6)
# ax[1].annotate(r'(a)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[1].legend(loc='upper right', fontsize=12)

# ax1 = ax[1].inset_axes([0.35, 0.2, 0.45, 0.45])

# linewidth_s = 1.
# markersize_s = 4
# fontsize_s = 16
# labelsize_s = 12


# diff = gf_imag_analytic - gf

# diff_imag = gf_imag_analytic - gf_imag

# ax1.plot(ts, diff, linewidth=linewidth_s, color=color1, markersize=markersize_s, markerfacecolor='none', ls='--')

# ax1.plot(ts_imag, diff_imag, linewidth=linewidth_s, color=color0, markersize=markersize_s, markerfacecolor='none', ls='--')


# ax1.set_ylabel(r'Error', fontsize=fontsize_s)
# ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
# ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax1.locator_params(axis='both', nbins=5)



plt.tight_layout(pad=0.5)

# plt.savefig('toulouse2.pdf', bbox_inches='tight')

plt.show()


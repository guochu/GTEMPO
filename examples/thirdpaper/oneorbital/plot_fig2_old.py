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


def read_real_tempo(beta, t, U, order=7):
	mu = U/2
	t2 = 100.
	lb = -2.
	ub = 2.
	dt = 0.05
	filename = 'result/anderson_tempo1_beta%s_t%s_U%s_mu%s_tep%s_lb%s_ub%s_dt0.05_e%s.json'%(beta, t+20., U, mu, t2, lb, ub, order)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)

	return data['ts'], data['gf'], data['ws'], data['Aw']

def read_real_tempo_2(beta, t, U, order=7):
	mu = U/2
	t2 = 100.
	lb = -2.
	ub = 2.
	dt = 0.05
	filename = 'result/anderson_tempo1_beta%s_t%s_U%s_mu%s_tep%s_lb%s_ub%s_dt0.05_e%s_2.json'%(beta, t+20., U, mu, t2, lb, ub, order)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)

	return data['ts'], data['gf'], data['ws'], data['Aw']

def read_imag_tempo(beta, U, dt=0.1, order=10):
	N = round(beta/dt)
	mu = U/2
	filename = 'result/anderson_tempo1_beta%s_U%s_mu%s_N%s_dt%s_imag_e%s.json'%(beta, U, mu, N, dt, order)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], -asarray(data['gf'])

def read_mc_data(filename):
	data = loadtxt(filename)
	return data[:, 0], data[:, 1]


fontsize = 20
labelsize = 16
linewidth = 1.5
markersize = 6

colors = ['b', 'c', 'y']
markers = ['o', '^', '+']


color0 = 'g'
color1 = 'b'
marker1 = 'o'
marker2 = '+'

fig, ax = plt.subplots(2, 2, figsize=(8, 7))



beta = 40.

linewidths = [2, 1.5, 1]
# alphas = [1, 0.7, 0.4]
alphas =  [0.3, 0.7, 1.]
ts = [5., 10., 15.]

# the first row with U=0.1
U = 0.1

real_order = 7

for i in range(3):
	t = ts[i]
	_linewidth = linewidths[i]
	alpha = alphas[i]


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo(beta, t, U, real_order)
	ts_e7_2, gf_e7_2, ws_e7_2, Aw_e7_2 = read_real_tempo_2(beta, t, U, real_order)

	ax[0,0].plot(ws_e7_2, Aw_e7_2, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))

	ax[0,1].plot(ts_e7_2, gf_e7_2, ls='--', color=color1, alpha=alpha, linewidth=linewidth)

ts_imag, gf_imag = read_imag_tempo(beta, U, 0.1)
ax[0,1].plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/paper03/cthyb/1orb/U%s/beta40/G-7.dat'%(U)
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax[0,1].plot(ts1, gf1, linewidth=0.5, color='k', ls='--', alpha=0.7, label=r'CTQMC')

# ax[0,0].legend(fontsize=12)
ax[0,0].set_ylabel(r'$A(\omega)$', fontsize=fontsize)
ax[0,0].set_xlabel(r'$\omega$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].annotate(r'(a)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(fontsize=12)

ax[0,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='x', nbins=6)
ax[0,1].annotate(r'(b)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)


ax1 = ax[0,1].inset_axes([0.3, 0.25, 0.4, 0.4])

linewidth_s = 1.
markersize_s = 4
fontsize_s = 16
labelsize_s = 12


for i in range(3):
	t = ts[i]
	_linewidth = linewidths[i]
	alpha = alphas[i]


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo(beta, t, U, real_order)
	ts_e7_2, gf_e7_2, ws_e7_2, Aw_e7_2 = read_real_tempo_2(beta, t, U, real_order)


	ax1.plot(ts_e7_2, gf_e7_2, ls='--', color=color1, alpha=alpha, linewidth=linewidth)


ax1.plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)
ax1.plot(ts1, gf1, linewidth=1, color='k', ls='--', alpha=0.7, label=r'CTQMC')

ax1.set_xlim(15, 25)
ax1.set_ylim(bottom=-0.205, top=-0.185)

# ax1.set_xlim(0, 0.5)
# ax1.set_ylim(bottom=-0.5, top=-0.46)


ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax1.locator_params(axis='both', nbins=3)



# the second row with U=0.5
U = 1.

for i in range(3):
	t = ts[i]
	_linewidth = linewidths[i]
	alpha = alphas[i]


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo(beta, t, U, real_order)
	ts_e7_2, gf_e7_2, ws_e7_2, Aw_e7_2 = read_real_tempo_2(beta, t, U, real_order)

	ax[1,0].plot(ws_e7_2, Aw_e7_2, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))

	ax[1,1].plot(ts_e7_2, gf_e7_2, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))

ts_imag, gf_imag = read_imag_tempo(beta, U, 0.1)
ax[1,1].plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/paper03/cthyb/1orb/U%s/beta40/G-7.dat'%(round(U))
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax[1,1].plot(ts1, gf1, linewidth=1, color='k', ls='--', alpha=0.7, label=r'CTQMC')


ax[1,0].set_ylabel(r'$A(\omega)$', fontsize=fontsize)
ax[1,0].set_xlabel(r'$\omega$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].annotate(r'(c)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)


ax[1,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='x', nbins=6)
ax[1,1].annotate(r'(d)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)


ax1 = ax[1,1].inset_axes([0.3, 0.25, 0.4, 0.4])

linewidth_s = 1.
markersize_s = 4
fontsize_s = 16
labelsize_s = 12


for i in range(3):
	t = ts[i]
	_linewidth = linewidths[i]
	alpha = alphas[i]


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo(beta, t, U, real_order)
	ts_e7_2, gf_e7_2, ws_e7_2, Aw_e7_2 = read_real_tempo_2(beta, t, U, real_order)

	ax1.plot(ts_e7_2, gf_e7_2, ls='--', color=color1, alpha=alpha, linewidth=linewidth)


ax1.plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)
ax1.plot(ts1, gf1, linewidth=1, color='k', ls='--', alpha=0.7, label=r'CTQMC')

ax1.set_xlim(15, 25)
ax1.set_ylim(bottom=-0.056, top=-0.049)
# ax1.set_ylabel(r'Error', fontsize=fontsize_s)
ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax1.locator_params(axis='both', nbins=5)


plt.tight_layout(pad=0.5)

plt.savefig('oneorbital2.pdf', bbox_inches='tight')

plt.show()

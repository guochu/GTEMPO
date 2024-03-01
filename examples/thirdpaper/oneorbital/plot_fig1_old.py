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

def read_real_tempo(beta, t0, U, dt, order):
	mu = U/2
	# t2 = t/2
	t = t0 + 20.
	filename = 'result/anderson_tempo1_beta%s_t%s_%s_U%s_mu%s_dt%s_e%s.json'%(beta, t, t0, U, mu, dt, order)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gf = -1j * gt + 1j * lt
	# lt = asarray(data['lt'])
	# gf_ts = dt * asarray(data['gf_ts'])
	gf_ts = dt * asarray([i for i in range(len(data['gf_ts']))])

	return data['ts'], asarray(data['ns']), gf_ts, gf, data['bd']

def read_real_tempo_2(beta, t0, U, dt, order):
	mu = U/2
	# t2 = t/2
	t = t0 + 20.
	filename = 'result/anderson_tempo1_beta%s_t%s_%s_U%s_mu%s_dt%s_e%s_2.json'%(beta, t, t0, U, mu, dt, order)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gf = -1j * gt + 1j * lt
	# lt = asarray(data['lt'])
	# gf_ts = dt * asarray(data['gf_ts'])
	gf_ts = dt * asarray([i for i in range(len(data['gf_ts']))])

	return data['ts'], asarray(data['ns']), gf_ts, gf, data['bd']


def read_real_tempo_b(beta, t, U, order=7):
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

def read_real_tempo_b_2(beta, t, U, order=7):
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
linewidth = 2.
markersize = 6

colors = ['b', 'c', 'g']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'b'

fig, ax = plt.subplots(3, 2, figsize=(7, 8))


beta = 40.
dt = 0.05

t = 20.
U = 0.5
real_order = 7

ts_e7, ns_e7, gf_ts_e7, gf_e7, bds_e7 = read_real_tempo(beta, t, U, dt, real_order)
ts_e7_2, ns_e7_2, gf_ts_e7_2, gf_e7_2, bds_e7_2 = read_real_tempo_2(beta, t, U, dt, real_order)

ax[0,0].plot(ts_e7, ns_e7, ls='-', color=color0, linewidth=linewidth, label=r'vacuum')
ax[0,0].plot(ts_e7_2, ns_e7_2, ls='--', color=color1, linewidth=linewidth, label=r'thermal')


ax[0,0].set_ylabel(r'$n$', fontsize=fontsize)
ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].annotate(r'(a)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(fontsize=12)

linewidths = [1.5, 1.5, 1.5]
# alphas = [1, 0.7, 0.3]
alphas =  [0.3, 0.7, 1.]
ts = [5., 10., 15.]

for i in range(3):
	t = ts[i]
	# _linewidth = linewidths[i]
	alpha = alphas[i]

	ts_e7, ns_e7, gf_ts_e7, gf_e7, bds_e7 = read_real_tempo(beta, t, U, dt, real_order)
	ts_e7_2, ns_e7_2, gf_ts_e7_2, gf_e7_2, bds_e7_2 = read_real_tempo_2(beta, t, U, dt, real_order)

	ax[0,1].plot(gf_ts_e7, -gf_e7.imag, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))
	# ax[0,1].plot(gf_ts_e7_2, -gf_e7_2.imag, ls='-', color=color1, alpha=alpha, linewidth=linewidth)

# print(abs(gf_e7.imag - gf_e7_2.imag).max())

ax[0,1].set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax[0,1].set_xlabel(r'$t-t_0$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='x', nbins=6)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].annotate(r'(b)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,1].legend(fontsize=12)

# vacuum state
for i in range(3):
	t = ts[i]
	_linewidth = linewidths[i]
	alpha = alphas[i]


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo_b(beta, t, U, real_order)

	ax[1,0].plot(ws_e7, Aw_e7, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))

	ax[1,1].plot(ts_e7, gf_e7, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))


ts_imag, gf_imag = read_imag_tempo(beta, U, 0.1)
ax[1,1].plot(ts_imag, gf_imag, ls='-', color='k', linewidth=0.5)

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/paper03/cthyb/1orb/U%s/beta40/G-7.dat'%(U)
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax[1,1].plot(ts1, gf1, linewidth=1, color='k', ls='--', alpha=0.7, label=r'CTQMC')

# ax[0,0].legend(fontsize=12)
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


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo_b(beta, t, U, real_order)

	ax1.plot(ts_e7, gf_e7, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))


ax1.plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)
ax1.plot(ts1, gf1, linewidth=1, color='k', ls='--', alpha=0.7, label=r'CTQMC')

# ax1.set_xlim(16, 24)
# ax1.set_ylim(bottom=-0.2, top=-0.185)

ax1.set_xlim(15, 25)
ax1.set_ylim(bottom=-0.13, top=-0.117)
# ax1.set_ylim(bottom=-0.21, top=-0.18)

ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax1.locator_params(axis='both', nbins=3)


# thermal state
for i in range(3):
	t = ts[i]
	_linewidth = linewidths[i]
	alpha = alphas[i]


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo_b_2(beta, t, U, real_order)

	ax[2,0].plot(ws_e7, Aw_e7, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))

	ax[2,1].plot(ts_e7, gf_e7, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))


ts_imag, gf_imag = read_imag_tempo(beta, U, 0.1)
ax[2,1].plot(ts_imag, gf_imag, ls='-', color='k', linewidth=0.5)

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/paper03/cthyb/1orb/U%s/beta40/G-7.dat'%(U)
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax[2,1].plot(ts1, gf1, linewidth=1, color='k', ls='--', alpha=0.7, label=r'CTQMC')

# ax[0,0].legend(fontsize=12)
ax[2,0].set_ylabel(r'$A(\omega)$', fontsize=fontsize)
ax[2,0].set_xlabel(r'$\omega$', fontsize=fontsize)
ax[2,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,0].locator_params(axis='both', nbins=6)
ax[2,0].annotate(r'(e)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)


ax[2,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[2,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[2,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,1].locator_params(axis='x', nbins=6)
ax[2,1].annotate(r'(f)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)


ax1 = ax[2,1].inset_axes([0.3, 0.25, 0.4, 0.4])

linewidth_s = 1.
markersize_s = 4
fontsize_s = 16
labelsize_s = 12


for i in range(3):
	t = ts[i]
	_linewidth = linewidths[i]
	alpha = alphas[i]


	ts_e7, gf_e7, ws_e7, Aw_e7 = read_real_tempo_b(beta, t, U, real_order)

	ax1.plot(ts_e7, gf_e7, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t)))


ax1.plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)
ax1.plot(ts1, gf1, linewidth=1, color='k', ls='--', alpha=0.7, label=r'CTQMC')

# ax1.set_xlim(16, 24)
# ax1.set_ylim(bottom=-0.2, top=-0.185)

ax1.set_xlim(15, 25)
ax1.set_ylim(bottom=-0.13, top=-0.117)
# ax1.set_ylim(bottom=-0.21, top=-0.18)

ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax1.locator_params(axis='both', nbins=3)

plt.tight_layout(pad=0.5)

plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


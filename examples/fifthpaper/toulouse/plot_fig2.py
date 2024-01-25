import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace
import math


rc('text', usetex=True)


def read_partial(mu, beta, dt, N, order):
	filename = 'result/partial_if_mu%s_beta%s_dt%s_N%s_order%s.json'%(mu, beta, dt, N, order)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# lt = asarray(data['lt'])
	return data['ts'], asarray(data['ns']), gt, data['bd'], data['time']

def read_ti(mu, beta, dt, N, order, k=5):
	filename = 'result/ti_if_mu%s_beta%s_dt%s_N%s_order%s_prony1.0e-5_k%s.json'%(mu, beta, dt, N, order, k)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# lt = asarray(data['lt'])
	return data['ts'], asarray(data['ns']), gt, data['bd'], data['time']


def read_real_analytic(mu, dt, N):
	filename = 'result/analytic_mu%s_dt%s_N%s.json'%(mu, dt, N)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], asarray(data['gf'])


# print(gf_real_analytic[-20:])
# print(gf[-20:])

# ts, gf, bds = read_imag_tempo(40., 0., 0.1, 8)

# print(gf)

fontsize = 20
labelsize = 16
linewidth = 2
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

ts = [10.,20.,30.,40.,50.,60.]
beta = 20.
mu = 0
dt = 0.05
Ns = [round(t / dt) for t in ts]
order = 7

# ts_real_analytic, gf_real_analytic = read_real_analytic(mu, dt, N)

r_pa = [read_partial(mu, beta, dt, N, order) for N in Ns]
ts_pa = [item[-1] for item in r_pa]

print('partial If times ', ts_pa)

r_ti = [read_ti(mu, beta, dt, N, order) for N in Ns]
ts_ti = [item[-1] for item in r_ti]

print('TI If takes ', ts_ti)

# ax.plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')

# ax.plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[0], markersize=markersize, linewidth=linewidth, label=r'Partial IF')
# ax.plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='--', color=colors[1], markersize=markersize, linewidth=linewidth, label=r'TI IF')

ax.plot(ts, ts_pa, ls='--', color=colors[0], markersize=markersize, linewidth=linewidth, label=r'Partial IF')
ax.plot(ts, ts_ti, ls='--', color=colors[1], markersize=markersize, linewidth=linewidth, label=r'TI IF')


# ts_real_tempo_e7, ns_real_tempo_e7, gf_real_tempo_e7, bds_real_tempo_e7 = read_real_tempo(beta, t, mu, dt, 7)
# ax[0,0].plot(ts_real_tempo_e7, gf_real_tempo_e7, ls='--', color=colors[1], linewidth=linewidth, label=r'$\varsigma=10^{-7}$')

# # ts_real_tempo_e8, ns_real_tempo_e8, gf_real_tempo_e8, bds_real_tempo_e8 = read_real_tempo(beta, t, mu, dt, 8)
# # ax[0,0].plot(ts_real_tempo_e8, gf_real_tempo_e8, ls='--', color=colors[2], linewidth=linewidth, label=r'$\varsigma=10^{-8}$')


ax.set_ylabel(r'run time', fontsize=fontsize)
ax.set_xlabel(r'$t$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.locator_params(axis='both', nbins=6)
# ax.annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax.legend(fontsize=12)




plt.tight_layout(pad=0.5)

# plt.savefig('toulouse1.pdf', bbox_inches='tight')

plt.show()





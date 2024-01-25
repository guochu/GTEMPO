import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt
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
linewidth = 1.5
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(1,2, figsize=(8,4))

t = 30.
beta = 20.
mu = 0
dt = 0.05
N = round(t / dt)

ts_real_analytic, gf_real_analytic = read_real_analytic(mu, dt, N)
# ax[0,0].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')


for i, order in enumerate([6,7,8]):
	ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_partial(mu, beta, dt, N, order)
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_ti(mu, beta, dt, N, order)

	# ax[0,0].plot(ts_pa, gf_pa, ls='--', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$\varsigma=10^{-%s}$'%(order))
	# ax[0,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth)


	ax[0].plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[i], markersize=markersize, linewidth=linewidth)
	ax[0].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth,label=r'$\varsigma=10^{-%s}$'%(order))

ax[0].set_ylabel(r'Error', fontsize=fontsize)
ax[0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
ax[0].annotate(r'(a)', xy=(0.15, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].legend(loc='upper right', fontsize=12)


ax[0].legend(fontsize=14)


order = 7

ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_partial(mu, beta, dt, N, order)
ax[1].plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[1], markersize=markersize, linewidth=linewidth)

for i, k in enumerate([5,10,15]):
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_ti(mu, beta, dt, N, order, k)

	# ax[1,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

	ax[1].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

ax[1].set_ylabel(r'Error', fontsize=fontsize)
ax[1].set_xlabel(r'$t$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].locator_params(axis='both', nbins=6)
ax[1].annotate(r'(b)', xy=(0.15, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1].legend(loc='upper right', fontsize=12)

ax[1].legend(fontsize=14)

# ax.plot(ts_pa, gf_pa, ls='--', color=colors[0], markersize=markersize, linewidth=linewidth, label=r'Partial IF')
# ax.plot(ts_ti, gf_ti, ls='--', color=colors[1], markersize=markersize, linewidth=linewidth, label=r'TI IF')


# ts_real_tempo_e7, ns_real_tempo_e7, gf_real_tempo_e7, bds_real_tempo_e7 = read_real_tempo(beta, t, mu, dt, 7)
# ax[0,0].plot(ts_real_tempo_e7, gf_real_tempo_e7, ls='--', color=colors[1], linewidth=linewidth, label=r'$\varsigma=10^{-7}$')

# # ts_real_tempo_e8, ns_real_tempo_e8, gf_real_tempo_e8, bds_real_tempo_e8 = read_real_tempo(beta, t, mu, dt, 8)
# # ax[0,0].plot(ts_real_tempo_e8, gf_real_tempo_e8, ls='--', color=colors[2], linewidth=linewidth, label=r'$\varsigma=10^{-8}$')



# ax[0,0].set_ylim(top=1.2)



plt.tight_layout(pad=0.5)

# plt.savefig('toulouse1.pdf', bbox_inches='tight')

plt.show()





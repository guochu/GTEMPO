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

def read_ti(mu, beta, dt, N, order, prony=5, k=5):
	filename = 'result/ti_if_mu%s_beta%s_dt%s_N%s_order%s_prony%s_k%s.json'%(mu, beta, dt, N, order, prony, k)
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


def read_lrz_partial(mu, beta, dt, N, order):
	filename = 'result/lrz_partial_if_mu%s_beta%s_dt%s_N%s_order%s.json'%(mu, beta, dt, N, order)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# lt = asarray(data['lt'])
	return data['ts'], asarray(data['ns']), gt, data['bd'], data['time']

def read_lrz_ti(mu, beta, dt, N, order, prony=5, k=5):
	filename = 'result/lrz_ti_if_mu%s_beta%s_dt%s_N%s_order%s_prony%s_k%s.json'%(mu, beta, dt, N, order, prony, k)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# lt = asarray(data['lt'])
	return data['ts'], asarray(data['ns']), gt, data['bd'], data['time']


def read_lrz_real_analytic(mu, dt, N):
	filename = 'result/lrz_analytic_mu%s_dt%s_N%s.json'%(mu, dt, N)
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

fig, ax = plt.subplots(2, 2, figsize=(8, 7))

t = 60.
beta = 20.
mu = 0
dt = 0.05
N = round(t / dt)
ts = [10.,20.,30.,40.,50.,60.]
Ns = [round(t / dt) for t in ts]
order = 7




ts_real_analytic, gf_real_analytic = read_real_analytic(mu, dt, N)
ax[0,0].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')

ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_partial(mu, beta, dt, N, order)
# print('partial If takes ', t_pa)

ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_ti(mu, beta, dt, N, order)
# print('TI If takes ', t_ti)

ax[0,0].plot(ts_pa, gf_pa, ls='--', color=colors[0], markersize=markersize, linewidth=linewidth, label=r'Partial')
ax[0,0].plot(ts_ti, gf_ti, ls='--', color=colors[1], markersize=markersize, linewidth=linewidth, label=r'TTI')

ax[0,0].set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].set_title(r'Semi-circle', fontsize=fontsize)
ax[0,0].annotate(r'(a)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)

linewidth_s = 1.
markersize_s = 4
fontsize_s = 16
labelsize_s = 12

ax1 = ax[0,0].inset_axes([0.45, 0.4, 0.5, 0.5])

ax1.plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[0], markersize=markersize, linewidth=linewidth_s, label=r'Partial')
ax1.plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='--', color=colors[1], markersize=markersize, linewidth=linewidth_s, label=r'TTI')

ax1.set_ylabel(r'Error', fontsize=fontsize_s)
ax1.set_xlabel(r'$t$', fontsize=fontsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax1.legend(fontsize=12)


r_pa = [read_partial(mu, beta, dt, N, order) for N in Ns]
ts_pa = [item[-1] for item in r_pa]

r_ti = [read_ti(mu, beta, dt, N, order) for N in Ns]
ts_ti = [item[-1] for item in r_ti]

ax[1,0].loglog(ts, ts_pa, ls='--', color=colors[0], marker=markers[0], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'Partial')
ax[1,0].loglog(ts, ts_ti, ls='--', color=colors[1], marker=markers[1], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'TTI')


ax[1,0].set_ylabel(r'Run time (s)', fontsize=fontsize)
ax[1,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[1,0].locator_params(axis='both', nbins=6)
# ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].annotate(r'(c)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)


# lorentzian spectrum
t = 60.
N = round(t / dt)

ts_real_analytic, gf_real_analytic = read_lrz_real_analytic(mu, dt, N)
ax[0,1].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')


ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_lrz_partial(mu, beta, dt, N, order)
# # # print('partial If takes ', t_pa)

ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_lrz_ti(mu, beta, dt, N, order)

# ax[0,1].plot(ts_pa, gf_pa, ls='--', color=colors[0], markersize=markersize, linewidth=linewidth, label=r'Partial IF')
ax[0,1].plot(ts_ti, gf_ti, ls='--', color=colors[1], markersize=markersize, linewidth=linewidth, label=r'TTI')


ax[0,1].set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax[0,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].set_title(r'Lorentzian', fontsize=fontsize)
ax[0,1].annotate(r'(b)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)


ax1 = ax[0,1].inset_axes([0.45, 0.4, 0.5, 0.5])

ax1.plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[0], markersize=markersize, linewidth=linewidth_s, label=r'Partial')
ax1.plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='--', color=colors[1], markersize=markersize, linewidth=linewidth_s, label=r'TTI')

ax1.set_ylabel(r'Error', fontsize=fontsize_s)
ax1.set_xlabel(r'$t$', fontsize=fontsize_s)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

ax1.legend(fontsize=12)


ts = [10.,20.,30.,40.,50.,60.]
Ns = [round(t / dt) for t in ts]

r_pa = [read_lrz_partial(mu, beta, dt, N, order) for N in Ns]
ts_pa = [item[-1] for item in r_pa]

ax[1,1].loglog(ts, ts_pa, ls='--', color=colors[0], marker=markers[0], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'Partial')


r_ti = [read_lrz_ti(mu, beta, dt, N, order) for N in Ns]
ts_ti = [item[-1] for item in r_ti]

ax[1,1].loglog(ts, ts_ti, ls='--', color=colors[1], marker=markers[1], markersize=markersize, markerfacecolor='none', linewidth=linewidth, label=r'TTI')


ax[1,1].set_ylabel(r'Run time (s)', fontsize=fontsize)
ax[1,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[1,1].locator_params(axis='both', nbins=6)
# ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].annotate(r'(d)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)



plt.tight_layout(pad=0.5)

plt.savefig('toulouse1.pdf', bbox_inches='tight')

plt.show()





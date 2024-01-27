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
linewidth = 1.5
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(3,2, figsize=(6,8))

t = 30.
beta = 20.
mu = 0
dt = 0.05
N = round(t / dt)

# semi-circle
ts_real_analytic, gf_real_analytic = read_real_analytic(mu, dt, N)
# ax[0,0].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')


for i, order in enumerate([6,7,8]):
	ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_partial(mu, beta, dt, N, order)
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_ti(mu, beta, dt, N, order)

	# ax[0,0].plot(ts_pa, gf_pa, ls='--', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$\varsigma=10^{-%s}$'%(order))
	# ax[0,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth)


	ax[2,0].plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[i], markersize=markersize, linewidth=linewidth)
	ax[2,0].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth,label=r'$\varsigma=10^{-%s}$'%(order))

ax[2,0].set_ylabel(r'Error', fontsize=fontsize)
ax[2,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[2,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,0].locator_params(axis='both', nbins=6)
ax[2,0].annotate(r'(e)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)
ax[2,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,0].legend(loc='upper right', fontsize=12)


# ax[0,0].legend(fontsize=14)


order = 7
prony = 5

ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_partial(mu, beta, dt, N, order)
ax[1,0].plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[1], markersize=markersize, linewidth=linewidth)

for i, k in enumerate([5,10,15]):
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_ti(mu, beta, dt, N, order, prony, k)

	# ax[1,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

	ax[1,0].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

ax[1,0].set_ylabel(r'Error', fontsize=fontsize)
ax[1,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].annotate(r'(c)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].legend(loc='upper right', fontsize=12)

ax[1,0].legend(fontsize=12)

for i, prony in enumerate([3,4,5]):
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_ti(mu, beta, dt, N, order, prony, 5)

	# ax[1,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

	ax[0,0].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$\varsigma_p=10^{-%s}$'%(prony))

ax[0,0].set_ylabel(r'Error', fontsize=fontsize)
ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].annotate(r'(a)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,0].legend(loc='upper right', fontsize=12)
ax[0,0].set_title(r'Semi-circle', fontsize=fontsize)

ax[0,0].legend(fontsize=12)


# Lorentzian

t = 30.
beta = 20.
mu = 0
dt = 0.05
N = round(t / dt)

# semi-circle
ts_real_analytic, gf_real_analytic = read_lrz_real_analytic(mu, dt, N)
# ax[0,0].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')


for i, order in enumerate([6,7,8]):
	ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_lrz_partial(mu, beta, dt, N, order)
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_lrz_ti(mu, beta, dt, N, order)

	# ax[0,0].plot(ts_pa, gf_pa, ls='--', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$\varsigma=10^{-%s}$'%(order))
	# ax[0,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth)


	ax[2,1].plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[i], markersize=markersize, linewidth=linewidth)
	ax[2,1].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth,label=r'$\varsigma=10^{-%s}$'%(order))

ax[2,1].set_ylabel(r'Error', fontsize=fontsize)
ax[2,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[2,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[2,1].locator_params(axis='both', nbins=6)
ax[2,1].annotate(r'(f)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)
ax[2,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[2,1].legend(loc='upper right', fontsize=12)


# ax[0,0].legend(fontsize=14)


order = 7
prony = 5

ts_pa, ns_pa, gf_pa, bd_pa, t_pa = read_lrz_partial(mu, beta, dt, N, order)
ax[1,0].plot(ts_pa, abs(gf_pa - gf_real_analytic[1:]), ls='--', color=colors[1], markersize=markersize, linewidth=linewidth)

for i, k in enumerate([5,10,15]):
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_lrz_ti(mu, beta, dt, N, order, prony, k)

	# ax[1,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

	ax[1,1].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

ax[1,1].set_ylabel(r'Error', fontsize=fontsize)
ax[1,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].annotate(r'(d)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)
ax[1,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,1].legend(loc='upper right', fontsize=12)

ax[1,1].legend(fontsize=12)

for i, prony in enumerate([3,4,5]):
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti = read_lrz_ti(mu, beta, dt, N, order, prony, 5)

	# ax[1,0].plot(ts_ti, gf_ti, ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$k=%s$'%(k))

	ax[0,1].plot(ts_ti, abs(gf_ti - gf_real_analytic[1:]), ls='-', color=colors[i], markersize=markersize, linewidth=linewidth, label=r'$\varsigma_p=10^{-%s}$'%(prony))

ax[0,1].set_ylabel(r'Error', fontsize=fontsize)
ax[0,1].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].annotate(r'(b)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)
ax[0,1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0,1].legend(loc='upper right', fontsize=12)
ax[0,1].set_title(r'Lorentzian', fontsize=fontsize)

ax[0,1].legend(fontsize=12)





plt.tight_layout(pad=0.5)

plt.savefig('toulouse2.pdf', bbox_inches='tight')

plt.show()





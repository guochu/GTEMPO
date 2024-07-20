import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, sqrt
from numpy.linalg import norm
import math


rc('text', usetex=True)


def read_partial(mu, beta, dt, N, order, chi=50):
	filename = 'result/partial_if_mu%s_beta%s_dt%s_N%s_order%s_chi%s.json'%(mu, beta, dt, N, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# lt = asarray(data['lt'])
	return data['ts'], asarray(data['ns']), gt, data['bd'], data['time'], asarray(data['S'])

def read_ti(mu, beta, dt, N, order, prony=5, k=5, chi=50):
	filename = 'result/ti_if_mu%s_beta%s_dt%s_N%s_order%s_prony%s_k%s_chi%s.json'%(mu, beta, dt, N, order, prony, k, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# lt = asarray(data['lt'])
	return data['ts'], asarray(data['ns']), gt, data['bd'], data['time'], asarray(data['S'])


def read_real_analytic(mu, dt, N):
	filename = 'result/analytic_mu%s_dt%s_N%s.json'%(mu, dt, N)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], asarray(data['gf'])

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

# def mse_error(a, b):
# 	assert len(a) == len(b)
# 	L = len(a)
# 	diff = asarray(a) - asarray(b)
# 	return abs(diff).sum() / L

fontsize = 20
labelsize = 16
linewidth = 1.5
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r', 'k']
markers = ['o', '^', '+']

fig, ax = plt.subplots(1,2, figsize=(8,4))


t = 60.
beta = 20.
mu = 0
dt = 0.05
N = round(t / dt)
ts = [10.,20.,30.,40.,50.,60.]
Ns = [round(t / dt) for t in ts]


# semi-circle
ts_real_analytic, gf_real_analytic = read_real_analytic(mu, dt, N)
# ax[0].plot(ts_real_analytic, gf_real_analytic, ls='-', color='k', linewidth=linewidth, label=r'Analytic')

# ts_pa, ns_pa, gf_pa, bd_pa, t_pa, S = read_partial(mu, beta, dt, N, order, chi)

# ax[0].plot(S, ls='-', color='k', linewidth=linewidth, label=r'Analytic')

# dt = 0.05
Sm = []
errs = []
dt = 0.05

order = 10
chi = 100

for j, tj in enumerate(ts):
	Nj = round(tj / dt)
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti, S = read_ti(mu, beta, dt, Nj, order, chi=chi)
	tj = round(ts[j])
	Sm.append(S.max())
	# ed_gf_adapted = gf_real_analytic[1:2*len(gf_ti):2]
	ed_gf_adapted = gf_real_analytic[1:len(gf_ti)+1]
	errs.append(mse_error(gf_ti, ed_gf_adapted))
	# ax[0].plot(linspace(0, 1, num=len(S)), S, ls='-', color=colors[j], linewidth=linewidth, label=r'$t=%s$'%(tj))

ax[0].plot(ts, Sm, ls='--', color='b', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none', label=r'$D\delta t=%s$'%(2*dt))

ax[0].set_ylabel(r'OSEE', fontsize=fontsize)
ax[0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
# ax[0].set_title(r'Lorentzian', fontsize=fontsize)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].set_title(r'$D\delta t=%s$'%(2*dt), fontsize=fontsize)
ax[0].annotate(r'(a)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)


ax1 = ax[0].inset_axes([0.4, 0.2, 0.5, 0.5])

ax1.plot(ts, errs, ls='--', color='k', linewidth=linewidth, markersize=markersize, markerfacecolor='none', label=r'$D\delta t=%s$'%(2*dt))

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax1.set_xlabel(r'$t$', fontsize=fontsize)
ax1.tick_params(axis='both', which='major', labelsize=labelsize)
ax1.locator_params(axis='both', nbins=6)
# ax[0].set_title(r'Lorentzian', fontsize=fontsize)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))


# dt = 0.1
Sm = []
errs = []
dt = 0.1

order = 10
chi = 100

for j, tj in enumerate(ts):
	Nj = round(tj / dt)
	ts_ti, ns_ti, gf_ti, bd_ti, t_ti, S = read_ti(mu, beta, dt, Nj, order, chi=chi)
	tj = round(ts[j])
	Sm.append(S.max())
	ed_gf_adapted = gf_real_analytic[1:2*len(gf_ti):2]
	errs.append(mse_error(gf_ti, ed_gf_adapted))
	# ax[0].plot(linspace(0, 1, num=len(S)), S, ls='-', color=colors[j], linewidth=linewidth, label=r'$t=%s$'%(tj))

ax[1].plot(ts, Sm, ls='--', color='b', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none', label=r'$D\delta t=%s$'%(2*dt))

ax[1].set_ylabel(r'OSEE', fontsize=fontsize)
ax[1].set_xlabel(r'$t$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].locator_params(axis='both', nbins=6)
# ax[0].set_title(r'Lorentzian', fontsize=fontsize)
ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1].set_title(r'$D\delta t=%s$'%(2*dt), fontsize=fontsize)
ax[1].annotate(r'(b)', xy=(-0.2, 1.05),xycoords='axes fraction', fontsize=fontsize)


ax1 = ax[1].inset_axes([0.4, 0.2, 0.5, 0.5])

ax1.plot(ts, errs, ls='--', color='k', linewidth=linewidth, markersize=markersize, markerfacecolor='none', label=r'$D\delta t=%s$'%(2*dt))

ax1.set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax1.set_xlabel(r'$t$', fontsize=fontsize)
ax1.tick_params(axis='both', which='major', labelsize=labelsize)
ax1.locator_params(axis='both', nbins=6)
# ax[0].set_title(r'Lorentzian', fontsize=fontsize)
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

# ax[0].legend(fontsize=12)


# ax[0].legend(fontsize=12)

plt.tight_layout(pad=0.5)

plt.savefig('toulouse3.pdf', bbox_inches='tight')

plt.show()

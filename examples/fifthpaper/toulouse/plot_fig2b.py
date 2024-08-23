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

def read_mixed_tempo(beta, dtau, t, mu, dt, order=10, chi=60):
	filename = 'result/thouless_tempo_mixed_beta%s_dtau%s_t%s_mu%s_dt%s_order%s_chi%s.json'%(beta, dtau, t, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gtau = parse_complex_array(data['gtau'])
	# gf =  gt + lt
	gf = 1j * gt + 1j * lt
	return data['taus'], -gtau.real, data['ts'], gt, lt, gf.imag, data['bd']

def read_ed(beta, N, mu, dt, dw=0.001):
	mpath = '../../../../iGTEMPO/examples/secondpaper/toulouse/'
	filename = mpath + 'result/Toulouse_ed_beta%s_mu%s_dt%s_N%s_dw%s.json'%(beta, mu, dt, N, dw)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	return data['ts'], gt, lt

def read_real_tempo(beta, t, mu, dt, order=10, chi=1024):
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/thouless_tempo_real_beta%s_t%s_mu%s_dt%s_order%s_chi%s.json'%(beta, t, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = asarray([item['re'] for item in data['gt']])
	# gt = parse_complex_array(data['gt'])
	# lt = parse_complex_array(data['lt'])
	# gf =  gt + lt
	return data['ts'], asarray(data['ns']), gt, data['bd']

def read_real_analytic(t, mu, dt):
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/thouless_analytic_real_t%s_mu%s_dt%s.json'%(t, mu, dt)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], asarray(data['gf'])


def read_imag_tempo(beta, mu, dt, order=10, chi=1024):
	N = round(beta/dt)
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/toulouse_tempo_beta%s_mu%s_N%s_order%s_chi%s.json'%(beta, mu, N, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	# lt = asarray(data['lt'])
	return data['ts'], -asarray(data['gf']), asarray(data['bd']).max()

def read_imag_analytic(beta, mu, dt):
	N = round(beta/dt)
	mpath = '../../thirdpaper/toulouse/'
	filename = mpath + 'result/thouless_analytic_beta%s_mu%s_N%s.json'%(beta, mu, N)
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
	return sqrt(v * v / L)
	# return abs(diff).max()

fontsize = 14
labelsize = 12
linewidth = 1.5
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2, 2, figsize=(8,6))

t = 5.
beta = 10.
mu = 0.
order = 12


chi = 120


dt = 0.05
all_dtau = [0.02, 0.1, 0.5]

errs_real = []
errs_imag = []


for (i, dtau) in enumerate(all_dtau):
	ts_real_analytic, gf_real_analytic = read_real_analytic(t, mu, dt=dt)
	ts_imag_analytic, gf_imag_analytic = read_imag_analytic(beta, mu, dt=dtau)

	taus, gtau, ts_tempo, gt_tempo, lt_tempo, gf_tempo, bds = read_mixed_tempo(beta, dtau, t, mu, dt, order=order, chi=chi)

	err_real = gf_real_analytic - gf_tempo
	errs_real.append(mse_error(gf_real_analytic, gf_tempo))

	ax[0,0].plot(ts_tempo, err_real, ls='--', color=colors[i], linewidth=linewidth, label=r'mixed, $\delta \tau=%s$'%(dtau))

	err_imag = gf_imag_analytic - gtau
	errs_imag.append(mse_error(gf_imag_analytic, gtau))

	ax[1,0].plot(taus, err_imag, ls='--', color=colors[i], linewidth=linewidth, label=r'mixed, $\delta \tau=%s$'%(dtau))


ax[0,0].set_ylabel(r'$-{\rm Im}[G_{Ana}(t) - G(t)]$', fontsize=fontsize)
ax[0,0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].annotate(r'(c)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(fontsize=12)


ax[1,0].set_ylabel(r'$G_{Ana}(\tau)-G(\tau)$', fontsize=fontsize)
ax[1,0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].annotate(r'(d)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[1,0].legend(fontsize=12)


ax[0,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[1,0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))



ax[0,1].loglog(all_dtau, errs_real, ls='--', color='k', linewidth=linewidth, markersize=markersize, marker='o', markerfacecolor='none')

ax[0,1].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[0,1].set_xlabel(r'$\delta \tau$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].annotate(r'(c)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)




ax[1,1].loglog(all_dtau, errs_imag, ls='--', color='k', linewidth=linewidth, markersize=markersize, marker='o', markerfacecolor='none')

ax[1,1].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
ax[1,1].set_xlabel(r'$\delta \tau$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[0,1].locator_params(axis='both', nbins=6)
ax[1,1].annotate(r'(d)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

# plt.savefig('toulouse2.pdf', bbox_inches='tight')

plt.show()

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


def read_real_tempo(beta, t, U, order=10, chi=60):
	mu = U/2
	t2 = 500.
	lb = -2.
	ub = 2.
	dt = 0.05
	filename = 'result/anderson_tempo1_beta%s_t%s_U%s_mu%s_tep%s_lb%s_ub%s_dt0.05_order%s_chi%s.json'%(beta, t+20., U, mu, t2, lb, ub, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)

	return data['ts'], data['gf'], data['ws'], data['Aw']

def read_real_tempo_2(beta, t, U, order=10, chi=60):
	mu = U/2
	t2 = 500.
	lb = -2.
	ub = 2.
	dt = 0.05
	filename = 'result/anderson_tempo1_beta%s_t%s_U%s_mu%s_tep%s_lb%s_ub%s_dt0.05_order%s_chi%s_2.json'%(beta, t+20., U, mu, t2, lb, ub, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)

	return data['ts'], data['gf'], data['ws'], data['Aw']

def read_imag_tempo(beta, U, dt=0.1, order=10, chi=500):
	N = round(beta/dt)
	mu = U/2
	filename = 'result/anderson_tempo1_beta%s_U%s_mu%s_N%s_dt%s_imag_order%s_chi%s.json'%(beta, U, mu, N, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], -asarray(data['gf'])

def read_mc_data(filename):
	data = loadtxt(filename)
	return data[:, 0], data[:, 1]

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)


fontsize = 20
labelsize = 16
linewidth = 1.5
markersize = 6

colors = ['b', 'c', 'y']
markers = ['o', '^', '+']


color0 = 'b'
color1 = 'c'
color2 = 'r'
marker1 = 'o'
marker2 = '+'

fig, ax = plt.subplots(3, 2, figsize=(6, 8))


beta = 40.

linewidths = [2, 1.5, 1]
# alphas = [1, 0.7, 0.4]
alphas =  [0.3, 0.7, 1.]
ts = [5., 15.]

# the first row with U=0.1
U = 0.1


t_small = 5.

times, gf, ws, Aw = read_real_tempo(beta, t_small, U)

alpha = 0.5

ax[0,0].plot(ws, Aw, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'Vacuum, $t_0=%s$'%(round(t_small)))
ax[0,1].plot(times, gf, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'Vacuum, $t_0=%s$'%(round(t_small)))

times, gf, ws, Aw = read_real_tempo_2(beta, t_small, U)

ax[0,0].plot(ws, Aw, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'Thermal, $t_0=%s$'%(round(t_small)))
ax[0,1].plot(times, gf, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'Thermal, $t_0=%s$'%(round(t_small)))


t_large = 80.

times, gf, ws, Aw = read_real_tempo_2(beta, t_large, U)

ax[0,0].plot(ws, Aw, ls='-', color=color2, alpha=alpha, linewidth=linewidth, label=r'Thermal, $t_0=%s$'%(round(t_large)))
ax[0,1].plot(times, gf, ls='-', color=color2, alpha=alpha, linewidth=linewidth, label=r'Thermal, $t_0=%s$'%(round(t_large)))




ts_imag, gf_imag = read_imag_tempo(beta, U, 0.1)
ax[0,1].plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1, label=r'imag, $\chi=500$')

mc_data_path = '/Users/guochu/Documents/Since2018/2024/MyPapers/gtempo/paper03/cthyb/1orb/U%s/beta40/G-7.dat'%(U)
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax[0,1].plot(ts1, gf1, linewidth=0.5, color='k', ls='--', alpha=0.7, label=r'CTQMC')


ax[0,0].set_ylabel(r'$A(\omega)$', fontsize=fontsize)
ax[0,0].set_xlabel(r'$\omega$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].locator_params(axis='both', nbins=6)
ax[0,0].annotate(r'(a)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,0].legend(fontsize=10)


ax[0,1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0,1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].locator_params(axis='x', nbins=6)
ax[0,1].annotate(r'(b)', xy=(0.05, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0,1].legend(fontsize=10)

# U = 0.5
U = 0.5

times, gf, ws, Aw = read_real_tempo(beta, t_small, U)

alpha = 0.5

ax[1,0].plot(ws, Aw, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t_small)))
ax[1,1].plot(times, gf, ls='--', color=color0, alpha=alpha, linewidth=linewidth)

times, gf, ws, Aw = read_real_tempo_2(beta, t_small, U)

ax[1,0].plot(ws, Aw, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t_small)))
ax[1,1].plot(times, gf, ls='--', color=color1, alpha=alpha, linewidth=linewidth)


times, gf, ws, Aw = read_real_tempo_2(beta, t_large, U)

ax[1,0].plot(ws, Aw, ls='-', color=color2, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t_large)))
ax[1,1].plot(times, gf, ls='-', color=color2, alpha=alpha, linewidth=linewidth)



ts_imag, gf_imag = read_imag_tempo(beta, U, 0.1)
ax[1,1].plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)

mc_data_path = '/Users/guochu/Documents/Since2018/2024/MyPapers/gtempo/paper03/cthyb/1orb/U%s/beta40/G-7.dat'%(U)
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax[1,1].plot(ts1, gf1, linewidth=0.5, color='k', ls='--', alpha=0.7, label=r'CTQMC')


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


# U = 1.
U = 1.

times, gf, ws, Aw = read_real_tempo(beta, t_small, U)

alpha = 0.5

ax[2,0].plot(ws, Aw, ls='--', color=color0, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t_small)))
ax[2,1].plot(times, gf, ls='--', color=color0, alpha=alpha, linewidth=linewidth)

times, gf, ws, Aw = read_real_tempo_2(beta, t_small, U)

ax[2,0].plot(ws, Aw, ls='--', color=color1, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t_small)))
ax[2,1].plot(times, gf, ls='--', color=color1, alpha=alpha, linewidth=linewidth)


times, gf, ws, Aw = read_real_tempo_2(beta, t_large, U)

ax[2,0].plot(ws, Aw, ls='-', color=color2, alpha=alpha, linewidth=linewidth, label=r'$t_0=%s$'%(round(t_large)))
ax[2,1].plot(times, gf, ls='-', color=color2, alpha=alpha, linewidth=linewidth)


ts_imag, gf_imag = read_imag_tempo(beta, U, 0.1)
ax[2,1].plot(ts_imag, gf_imag, ls='-', color='k', linewidth=1)

mc_data_path = '/Users/guochu/Documents/Since2018/2024/MyPapers/gtempo/paper03/cthyb/1orb/U%s/beta40/G-7.dat'%(round(U))
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax[2,1].plot(ts1, gf1, linewidth=0.5, color='k', ls='--', alpha=0.7, label=r'CTQMC')


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


plt.tight_layout(pad=0.5)

plt.savefig('oneorbital3.pdf', bbox_inches='tight')

plt.show()

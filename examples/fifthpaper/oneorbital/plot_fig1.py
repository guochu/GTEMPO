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

def read_mixed_tempo(beta, t, U, dt=0.05, order=10, chi=60):
	# mpath = '/Users/guochu/Documents/Missile/GTEMPO/examples/thirdpaper/oneorbital'
	mu = U/2
	# t2 = t/2
	filename = 'result/mixed_beta%s_dtau0.1_t%s_dt%s_U%s_mu%s_order%s_chi%s.json'%(beta, t, dt, U, mu, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gf = -1j * gt - 1j * lt
	ts = asarray(data['ts'])
	return ts-ts[0], gf, gt, lt


def read_real_tempo(beta, t0, U, dt, order=10, chi=60):
	mpath = '../../thirdpaper/oneorbital/'
	mu = U/2
	# t2 = t/2
	t = t0 + 20.
	filename = mpath + 'result/anderson_tempo1_beta%s_t%s_%s_U%s_mu%s_dt%s_order%s_chi%s_2.json'%(beta, t, t0, U, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	gf = -1j * gt + 1j * lt

	gf_ts = dt * asarray([i for i in range(len(data['gf_ts']))])

	return data['ts'], asarray(data['ns']), gf_ts, gf, gt, lt


def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)


fontsize = 20
labelsize = 16
linewidth = 2.
markersize = 10

colors = ['b', 'c', 'g']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'r'

fig, ax = plt.subplots(1, 1, figsize=(6, 5))


beta = 40.
dt = 0.05

# U = 0.1
U = 0.5

t_final = 40.

t0 = 20.

chi = 60

times_final, ns_final, gf_ts_final, gf_final, gt_final, lt_final = read_real_tempo(beta, t_final, U, dt, chi=chi)

ax.plot(gf_ts_final, -gf_final.imag, ls='--', color=color0, linewidth=linewidth, label=r'real, $t_0=%s, \chi=%s$'%(round(t_final), chi))

# print(lt_final[-10:])

chi = 100
mixed_ts, mixed_gf, mixed_gt, mixed_lt = read_mixed_tempo(beta, t0, U, dt, chi=chi)

# print(gf_ts_final)
# print(-mixed_lt[-10:])

ax.plot(mixed_ts, -mixed_gf.imag, ls='--', color=color1, linewidth=linewidth, label=r'mixed, $\chi=%s$'%(chi))

ax.set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax.set_ylabel(r'$-{\rm Im}[G^R(t)]$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')
ax.locator_params(axis='both', nbins=6)
# ax.set_ylim(0, 1)
# ax.set_xlim(0, t * 0.1)

ax.legend(fontsize=12)


plt.tight_layout(pad=0.5)

# plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()


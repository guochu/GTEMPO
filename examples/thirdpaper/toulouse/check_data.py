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

def read_real_tempo(beta, t, mu, dt, order=10, chi=1024):
	filename = 'result/thouless_tempo_real_beta%s_t%s_mu%s_dt%s_order%s_chi%s.json'%(beta, t, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	# lt = asarray(data['lt'])
	return data['ts'], gt, lt


def read_ed(beta, N, mu, dt, dw=0.001):
	mpath = '../../../../iGTEMPO/examples/secondpaper/toulouse/'
	filename = mpath + 'result/Toulouse_ed_beta%s_mu%s_dt%s_N%s_dw%s.json'%(beta, mu, dt, N, dw)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	return data['ts'], gt, lt

def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)


# print(gf_real_analytic[-20:])
# print(gf[-20:])

# ts, gf, bds = read_imag_tempo(40., 0., 0.1, 8)

# print(gf)

fontsize = 20
labelsize = 16
linewidth = 1.5
markersize = 7

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(2, 2, figsize=(8, 7))

t = 40.
beta = 40.
mu = 0.
dt = 0.05
N = round(t / dt)

ts_ed, gt_ed, lt_ed = read_ed(beta, N, mu, dt, dw=0.0005)
ax[0,0].plot(ts_ed, gt_ed.real, ls='--', color='k', linewidth=linewidth, label=r'ED')
ax[0,1].plot(ts_ed, gt_ed.imag, ls='--', color='k', linewidth=linewidth, label=r'ED')

ax[1,0].plot(ts_ed, lt_ed.real, ls='--', color='k', linewidth=linewidth, label=r'ED')
ax[1,1].plot(ts_ed, lt_ed.imag, ls='--', color='k', linewidth=linewidth, label=r'ED')


chi = 60

ts_tempo, gt_tempo, lt_tempo = read_real_tempo(beta, t, mu, dt, chi=chi)

ax[0,0].plot(ts_tempo, gt_tempo.real, ls='--', color='k', linewidth=linewidth, label=r'ED')
ax[0,1].plot(ts_tempo, gt_tempo.imag, ls='--', color='k', linewidth=linewidth, label=r'ED')

ax[1,0].plot(ts_tempo, lt_tempo.real, ls='--', color='k', linewidth=linewidth, label=r'ED')
ax[1,1].plot(ts_tempo, lt_tempo.imag, ls='--', color='k', linewidth=linewidth, label=r'ED')


plt.tight_layout(pad=0.5)

# plt.savefig('toulouse1.pdf', bbox_inches='tight')

plt.show()





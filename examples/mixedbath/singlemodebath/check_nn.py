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

def read_real_tempo(beta, t, N, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dt = t / N
	filename = 'result/noninteracting_realgtempo_beta%s_t%s_dt%s_omega0%s_alpha0%s_omega1%s_alpha1%s_mu%s_chi%s.json'%(beta, t, dt, omega0, alpha0, omega1, alpha1, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	nn2 = parse_complex_array(data['nn2'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt, nn, nn2

def read_mixed_tempo(beta, Ntau, t, N, mu, omega0, alpha0, omega1, alpha1, chi=80):
	dt = t / N
	dtau = beta / Ntau
	filename = 'result/noninteracting_mixedgtempo_beta%s_dtau%s_t%s_dt%s_omega0%s_alpha0%s_omega1%s_alpha1%s_mu%s_chi%s.json'%(beta, dtau, t, dt, omega0, alpha0, omega1, alpha1, mu, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	ts = asarray(data['ts'])
	return ts-ts[0], gt, lt, nn

def read_neq_ed(beta, t, N, mu, omega0, alpha0, omega1, alpha1):
	filename = 'result/noninteracting_neq_ED_real_beta%s_mu%s_t%s_N%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, t, N, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	return data['ts'], gt, lt, nn

def read_eq_ed(beta, t, N, mu, omega0, alpha0, omega1, alpha1):
	filename = 'result/noninteracting_eq_ED_real_beta%s_mu%s_t%s_N%s_omega0%s_alpha0%s_omega1%s_alpha1%s.json'%(beta, mu, t, N, omega0, alpha0, omega1, alpha1)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = parse_complex_array(data['gt'])
	lt = parse_complex_array(data['lt'])
	nn = parse_complex_array(data['nn'])
	return data['ts'], gt, lt, nn


def mse_error(a, b):
	assert len(a) == len(b)
	L = len(a)
	diff = asarray(a) - asarray(b)
	v = norm(diff)
	return sqrt(v * v / L)

fontsize = 18
labelsize = 14
linewidth = 2
markersize = 10

colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+']

fig, ax = plt.subplots(1,2, figsize=(8,4))


omega0 = 1
alpha0 = 0.25
omega1 = 1
alpha1 = 1
chi = 200

mu = 0.

t = 5
Nt = 400
beta = 5
# Ntau = 20

# real time data
ts_ed, gt_ed, lt_ed, nn_ed = read_neq_ed(beta, t, Nt, mu, omega0, alpha0, omega1, alpha1)
ts2, gt2, lt2, nn, nn2 = read_real_tempo(beta, t, Nt, mu, omega0, alpha0, omega1, alpha1, chi)

print(len(nn), ' ', len(nn2))
# print(nn[:5])
# print(nn2[:5])

ax[0].plot(ts_ed[2:], nn_ed.real[2:], ls='-', color='k', linewidth=linewidth, label=r'ED')
ax[0].plot(ts2[2:], nn.real, ls='--', color='r', linewidth=linewidth, markersize=markersize, markerfacecolor='none', label=r'new, $\chi=%s$'%(chi))
ax[0].plot(ts2[2:], nn2.real, ls='--', color='g', linewidth=linewidth, markersize=markersize, markerfacecolor='none', label=r'old, $\chi=%s$'%(chi))


print(nn.real[2:][:50])
print(nn2.real[:50])

ax[0].set_xlabel(r'$t$', fontsize=fontsize)
ax[0].set_ylabel(r'$nn$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].locator_params(axis='both', nbins=6)
ax[0].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax[0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].legend(loc='lower right', fontsize=10)


# chis = [100, 200,300, 400, 500, 600,700]
# nn_errs = []
# nn2_errs = []

# for chi in chis:
# 	ts2, gt2, lt2, nn, nn2 = read_real_tempo(beta, t, Nt, mu, omega0, alpha0, omega1, alpha1, chi)
# 	nn_errs.append(mse_error(nn_ed.real[2:], nn.real))
# 	nn2_errs.append(mse_error(nn_ed.real[2:], nn2.real))

# ax[1].plot(chis, nn_errs, ls='--', color='r', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none', label=r'new')
# ax[1].plot(chis, nn2_errs, ls='--', color='g', linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none', label=r'old')
# ax[1].set_xlabel(r'$\chi$', fontsize=fontsize)
# ax[1].set_ylabel(r'$\mathcal{E}$', fontsize=fontsize)
# ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
# ax[1].locator_params(axis='both', nbins=6)
# ax[1].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax[1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
# ax[1].legend(fontsize=10)

plt.tight_layout(pad=0.5)

# plt.savefig('independentbosons3.pdf', bbox_inches='tight')

plt.show()

import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace
import math, os


rc('text', usetex=True)


def read_tempo_data_e8(beta):
	N = int(beta / 0.1)
	chi=2048
	data_path = 'result/checkbond1_beta%s_N%s_chi%s.json'%(beta, N, chi)
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return asarray(data['bd'])


def read_tempo_data_e9(beta):
	N = int(beta / 0.1)
	chi=2048
	data_path = 'result/checkbond2_beta%s_N%s_chi%s.json'%(beta, N, chi)
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return asarray(data['bd'])

def read_toulouse_data(beta):
	mu = 1.
	N = int(beta / 0.1)
	chi=512
	data_path = '../thouless/result/thouless_tempo_beta%s_mu%s_N%s_chi%s.json'%(beta, mu, N, chi)
	if not os.path.isfile(data_path):
		chi=1024
		data_path = '../thouless/result/thouless_tempo_beta%s_mu%s_N%s_chi%s.json'%(beta, mu, N, chi)
	with open(data_path, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return asarray(data['bd'])

# betas_toulouse = [5., 10., 15., 20., 25., 30.]
betas_toulouse = [5., 10., 15., 20., 25., 30., 40., 50., 60.]
bds_toulouse = [read_toulouse_data(beta) for beta in betas_toulouse]
bds_toulouse = [item.max() for item in bds_toulouse]
print(bds_toulouse)


betas_e8 = [5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60.]
bds_e8 = [read_tempo_data_e8(beta) for beta in betas_e8]
bds_e8 = [item.max() for item in bds_e8]
print(bds_e8)

betas_e9 = [5., 10., 15., 20.]
bds_e9 = [read_tempo_data_e9(beta) for beta in betas_e9]
bds_e9 = [item.max() for item in bds_e9]
print(bds_e9)

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

fontsize = 18
labelsize = 14
linewidth = 1.
markersize = 10

ax.semilogy(betas_toulouse, [item*item for item in bds_toulouse], linewidth=linewidth, marker='^', markersize=markersize, markerfacecolor='none', color='g', ls='--', label=r'$\mathcal{I}_{\uparrow}\mathcal{I}_{\downarrow}, \varsigma=10^{-8}$')

ax.semilogy(betas_e9, bds_e9, linewidth=linewidth, marker='o', markersize=markersize, markerfacecolor='none', color='r', ls='--', label=r'$\tilde{\mathcal{I}}, \varsigma=10^{-9}$')

ax.semilogy(betas_e8, bds_e8, linewidth=linewidth, marker='x', markersize=markersize, markerfacecolor='none', color='b', ls='--', label=r'$\tilde{\mathcal{I}}, \varsigma=10^{-8}$')


ax.set_ylabel(r'$\chi$', fontsize=fontsize)
ax.set_xlabel(r'$\beta$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')

ax.legend(fontsize=labelsize)

plt.tight_layout(pad=0.5)

plt.savefig('oneorbital3.pdf', bbox_inches='tight')

plt.show()










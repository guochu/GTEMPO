import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace
import math


rc('text', usetex=True)

def read_tempo_data(filename):
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], -asarray(data['gf']), data['bd']

def read_mc_data(filename):
	data = loadtxt(filename)
	return data[:, 0], data[:, 1]


fontsize = 18
labelsize = 14
linewidth = 2.
markersize = 4

fig, ax = plt.subplots(1, 2, figsize=(8, 4))

U = 1.
mu = U/2
beta = 10.
N = 100

chi = 512
data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s.json'%(beta, U, mu, N, chi)
ts, gf, bds = read_tempo_data(data_path)

ax[0].plot(ts, gf, linewidth=linewidth, color='r', ls='-.', label=r'$\mathcal{K}\mathcal{I}_{\uparrow} \mathcal{I}_{\downarrow}, \varsigma=10^{-8}$')

print(asarray(bds).max())

chi = 1024

data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s_2.json'%(beta, U, mu, N, chi)
ts2, gf2, bds2 = read_tempo_data(data_path)

print(asarray(bds2).max())

ax[0].plot(ts2, gf2, linewidth=linewidth, color='y', ls='-.', label=r'$\mathcal{K}\tilde{\mathcal{I}}, \varsigma=10^{-8}$')


chi = 1024

data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s_2_e9.json'%(beta, U, mu, N, chi)
ts2, gf2, bds2 = read_tempo_data(data_path)

print(asarray(bds2).max())

ax[0].plot(ts2, gf2, linewidth=linewidth, color='b', ls='-.', label=r'$\mathcal{K}\tilde{\mathcal{I}}, \varsigma=10^{-9}$')

diff = abs((gf - gf2) / gf)
# print(diff)

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/chen/im/1orb/beta10/G-7.dat'
ts, gf = read_mc_data(mc_data_path)
ts = linspace(0, beta, num=len(ts))
ax[0].scatter(ts, gf, linewidth=0.2, marker='+', color='g', ls='--', label=r'CTMQC')



ax[0].legend(fontsize=12)

ax[0].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[0].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].tick_params(axis='y', which='both')
ax[0].annotate(r'(a)', xy=(-0.1, 1.05),xycoords='axes fraction', fontsize=fontsize)


beta = 20.
N = 200

chi = 512
data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s.json'%(beta, U, mu, N, chi)
ts, gf, bds = read_tempo_data(data_path)

ax[1].plot(ts, gf, linewidth=linewidth, color='r', ls='-.', label=r'$\mathcal{K}\mathcal{I}_{\uparrow} \mathcal{I}_{\downarrow}, \varsigma=10^{-8}$')

chi = 2048

data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s_2_e9.json'%(beta, U, mu, N, chi)
ts2, gf2, bds2 = read_tempo_data(data_path)

print(asarray(bds2).max())

ax[1].plot(ts2, gf2, linewidth=linewidth, color='b', ls='-.', label=r'$\mathcal{K}\tilde{\mathcal{I}}, \varsigma=10^{-9}$')

diff = abs((gf - gf2) / gf)
# print(diff)

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/chen/im/1orb/beta20/G-7.dat'
ts, gf = read_mc_data(mc_data_path)
ts = linspace(0, beta, num=len(ts))
ax[1].scatter(ts, gf, linewidth=0.2, marker='+', color='g', ls='--', label=r'CTMQC')



ax[1].legend(fontsize=12)

ax[1].set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax[1].set_xlabel(r'$\tau$', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].tick_params(axis='y', which='both')
ax[1].annotate(r'(b)', xy=(-0.1, 1.05),xycoords='axes fraction', fontsize=fontsize)


plt.tight_layout(pad=0.5)

plt.savefig('oneorbital2.pdf', bbox_inches='tight')

plt.show()

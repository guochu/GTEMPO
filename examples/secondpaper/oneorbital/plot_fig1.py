import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace
import math


rc('text', usetex=True)

def read_data(filename):
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], -asarray(data['gf'])

def read_mc_data(filename):
	data = loadtxt(filename)
	return data[:, 0], data[:, 1]


fontsize = 18
labelsize = 14
linewidth = 2.
markersize = 4

fig, ax = plt.subplots(1, 1, figsize=(6, 5))

U = 1.
mu = U/2

chi = 512

beta = 20.
N = 200
data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s.json'%(beta, U, mu, N, chi)
ts_tempo1, gf_tempo1 = read_data(data_path)

ax.plot([t/beta for t in ts_tempo1], gf_tempo1, linewidth=linewidth, color='g', ls='-.', label=r'$\beta=20$')


mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/chen/im/1orb/beta20/G-7.dat'
ts1, gf1 = read_mc_data(mc_data_path)
ts1 = linspace(0, beta, num=len(ts1))
ax.plot([t/beta for t in ts1], gf1, linewidth=1, color='g', ls='-')


beta = 40.
N = 400
data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s.json'%(beta, U, mu, N, chi)
ts_tempo2, gf_tempo2 = read_data(data_path)

ax.plot([t/beta for t in ts_tempo2], gf_tempo2, linewidth=linewidth, color='r', ls='-.', label=r'$\beta=40$')

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/chen/im/1orb/beta40/G-7.dat'
ts2, gf2 = read_mc_data(mc_data_path)
ts2 = linspace(0, beta, num=len(ts2))
ax.plot([t/beta for t in ts2], gf2, linewidth=1, color='r', ls='--')

beta = 60.
N = 600
data_path = 'result/anderson_tempo1_beta%s_U%s_e%s_N%s_chi%s.json'%(beta, U, mu, N, chi)
ts_tempo3, gf_tempo3 = read_data(data_path)

ax.plot([t/beta for t in ts_tempo3], gf_tempo3, linewidth=linewidth, color='y', ls='-.', label=r'$\beta=60$')

mc_data_path = '/Users/guochu/Documents/Since2018/2023/MyPapers/gtempo/chen/im/1orb/beta60/G-7.dat'
ts3, gf3 = read_mc_data(mc_data_path)
ts3 = linspace(0, beta, num=len(ts3))
ax.plot([t/beta for t in ts3], gf3, linewidth=1, color='y', ls='--')


ax.legend(fontsize=labelsize)

ax.set_ylabel(r'$G(\tau)$', fontsize=fontsize)
ax.set_xlabel(r'$\tau/\beta$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.tick_params(axis='y', which='both')


# ax1 = ax.inset_axes([0.3, 0.2, 0.5, 0.5])

# linewidth_s = 1.
# markersize_s = 4
# fontsize_s = 12
# labelsize_s = 10

# # diff1 = abs(gf - gf_tempo)

# ax1.plot([t/beta for t in ts1], diff1, linewidth=linewidth_s, color='g', alpha=0.2, markersize=markersize_s, markerfacecolor='none', ls='--', label=r'$\beta=20$')


# ax1.plot(ts, diff_e7, linewidth=linewidth_s, color='r', alpha=0.5, markersize=markersize_s, markerfacecolor='none', ls='--', label=r'GTEMPO')


# ax1.plot(ts, diff_e8, linewidth=linewidth_s, color='y', alpha=1, markersize=markersize_s, markerfacecolor='none', ls='--', label=r'GTEMPO')


# ax1.set_ylabel(r'Error', fontsize=fontsize_s)
# ax1.set_xlabel(r'$\tau$', fontsize=fontsize_s, labelpad=0.1)
# ax1.tick_params(axis='both', which='major', labelsize=labelsize_s)
# ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# ax1.locator_params(axis='both', nbins=6)



plt.tight_layout(pad=0.5)

plt.savefig('oneorbital1.pdf', bbox_inches='tight')

plt.show()

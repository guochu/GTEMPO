import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt
import math

import h5py

rc('text', usetex=True)

def rescale(curr):
	return asarray([2*math.pi*item/0.1 for item in curr])

def read_ed_data_single(filename):
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], data['Ileft'], data['Iright']

def read_ed_data(Vs, dw):
	Ilefts = []
	Irights = []
	for V in Vs:
		filename = 'result/anderson_ed_V%s_dw%s.json'%(V, dw)
		ts, Ileft, Iright = read_ed_data_single(filename)
		Ilefts.append(Ileft[-1])
		Irights.append(Iright[-1])
	return rescale(Ilefts), rescale(Irights)

def read_tempo_data_single(filename):
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], data['Ileft'], data['Iright'], data['bd']


def read_tempo_data(Vs, U):
	Ilefts = []
	Irights = []
	bds = []
	for V in Vs:
		filename = 'result/anderson_tempo1_V%s_U%s_N600_b.json'%(V, U)
		ts, Ileft, Iright, bd = read_tempo_data_single(filename)
		Ilefts.append(Ileft[-1])
		Irights.append(Iright[-1])
		bds.append( asarray(bd).max() )
	return rescale(Ilefts), rescale(Irights), bds


def read_tempo_data2(Vs, U, dt, order=7, k=5, chi=160):
	t = 6.3
	N = round(t / dt)
	Ilefts = []
	Irights = []
	bds = []
	for V in Vs:
		filename = 'result/ti_V%s_U%s_N%s_dt%s_order%s_chi%s_k%s.json'%(V, U, N, dt, order, chi, k)
		ts, Ileft, Iright, bd = read_tempo_data_single(filename)
		Ilefts.append(rescale(Ileft))
		Irights.append(rescale(Iright))
		bds.append( asarray(bd).max() )
	return Ilefts, Irights, bds


Vs = [0.17857143, 0.35714286, 0.53571429, 0.71428571, 0.89285715]


dt = 0.014
order = 7
k = 5
chi = 160
# our data
U = 0.
Ilefts_0, Irights_0, bd_0 = read_tempo_data2(Vs, U, dt=dt, order=order, k=k, chi=chi)
Imean_0 = [0.5 * ( a - b ) for a, b in zip(Ilefts_0, Irights_0)]

print(Ilefts_0[4])
print(Irights_0[4])

# print([item[-1] for item in Ilefts_0])
# print([item[-1] for item in Irights_0])

# print(Imean_0[2])
print([item[0] for item in Imean_0])
print([item[-1] for item in Imean_0])

# U = 2.
# Ilefts_2, Irights_2, bd_2 = read_tempo_data(Vs2, U)
# Imean_2 = 0.5 * ( Ilefts_2 - Irights_2 )

# U = 4.
# Ilefts_4, Irights_4, bd_4 = read_tempo_data(Vs2, U)
# Imean_4 = 0.5 * ( Ilefts_4 - Irights_4 )


# U = 6.
# Ilefts_6, Irights_6, bd_6 = read_tempo_data(Vs2, U)
# Imean_6 = 0.5 * ( Ilefts_6 - Irights_6 )

# U = 8.
# Ilefts_8, Irights_8, bd_8 = read_tempo_data(Vs2, U)
# Imean_8 = 0.5 * ( Ilefts_8 - Irights_8 )


# print(Ilefts_8)
# print(Irights_8)
# print(asarray(Ilefts_8) + asarray(Irights_8))

# print(bd_2)
# print(bd_4)
# print(bd_6)
# print(bd_8)



# fontsize = 16
# labelsize = 14
# linewidth = 1.
# markersize = 6

# fig, ax = plt.subplots(1, 1, figsize=(6, 4.5))


# # U = 0.
# alpha = 1.

# ax.plot(Vs2, Imean_0, alpha=alpha, linewidth=1.5, color='b', marker='+', markersize=markersize, markerfacecolor='none', ls = 'none', label=r'GTEMPO')

# alpha = 0.8
# ax.plot(Vs2, Imean_2, linewidth=1.5, alpha=alpha, color='b', marker='+', markersize=markersize, markerfacecolor='none', ls = 'none')

# alpha = 0.6
# ax.plot(Vs2, Imean_4, linewidth=1.5, alpha=alpha, color='b', marker='+', markersize=markersize, markerfacecolor='none', ls = 'none')

# alpha = 0.4
# ax.plot(Vs2, Imean_6, linewidth=1.5, alpha=alpha, color='b', marker='+', markersize=markersize, markerfacecolor='none', ls = 'none')

# alpha = 0.2
# ax.plot(Vs2, Imean_8, linewidth=1.5, alpha=alpha, color='b', marker='+', markersize=markersize, markerfacecolor='none', ls = 'none')

# ax.text(2.,1.7,r'$U/\Gamma = 0$',rotation=30, fontsize=10)
# ax.text(2.,1.48,r'$U/\Gamma = 2$',rotation=30, fontsize=10)
# ax.text(2.,1.08,r'$U/\Gamma = 4$',rotation=25, fontsize=10)
# ax.text(2.,0.75,r'$U/\Gamma = 6$',rotation=20, fontsize=10)
# ax.text(2.,0.54,r'$U/\Gamma = 8$',rotation=15, fontsize=10)


# ax.set_ylabel(r'$2\pi \mathcal{J}/\Gamma$', fontsize=fontsize)
# ax.set_xlabel(r'$V/\Gamma$', fontsize=fontsize)
# ax.tick_params(axis='both', which='major', labelsize=labelsize)
# ax.tick_params(axis='y', which='both')
# ax.locator_params(axis='both', nbins=6)
# ax.set_ylim(0, 2.)
# ax.set_xlim(0, 3)
# ax.legend(fontsize=12)

# plt.tight_layout(pad=0.5)

# # plt.savefig('fig3.pdf', bbox_inches='tight')

# plt.show()

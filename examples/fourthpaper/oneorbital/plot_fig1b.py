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

def read_ed_data(Vs, dw=0.001):
	Ilefts = []
	Irights = []
	mpath = '/Users/guochu/Documents/QuantumSimulator/MPSImpuritySolvers/examples/firstpaper/twobath/'
	for V in Vs:
		filename = mpath + 'result/anderson_ed_V%s_dw%s.json'%(V, dw)
		ts, Ileft, Iright = read_ed_data_single(filename)
		Ilefts.append(Ileft[-1])
		Irights.append(Iright[-1])
	return rescale(Ilefts), rescale(Irights)


def read_ed_data_n(Vs, dw=0.001):
	Ilefts = []
	Irights = []
	for V in Vs:
		filename = 'result/anderson_ed_V%s_t6.3_dw%s.json'%(V, dw)
		ts, Ileft, Iright = read_ed_data_single(filename)
		Ilefts.append(Ileft[-1])
		Irights.append(Iright[-1])
	return rescale(Ilefts), rescale(Irights)

def read_tempo_data_single(filename):
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return data['ts'], data['Ileft'], data['Iright'], data['bd']


def read_tempo_data(Vs, U, dt=0.014, order=7):
	Ilefts = []
	Irights = []
	bds = []
	t = 4.2
	N = round(t/dt)
	mpath = '/Users/guochu/Documents/QuantumSimulator/MPSImpuritySolvers/examples/firstpaper/twobath/'
	for V in Vs:
		filename = mpath + 'result/anderson_tempo1_V%s_U%s_N%s_dt%s_order%s.json'%(V, U, N, dt, order)
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
colors = ['b', 'g', 'c', 'y', 'r']
markers = ['o', '^', '+', 'x', '*']

ed_Ilefts_old, ed_Irights_old = read_ed_data(Vs)
ed_Imean_old = 0.5 * ( ed_Ilefts_old - ed_Irights_old )

ed_Ilefts, ed_Irights = read_ed_data_n(Vs)
ed_Imean = 0.5 * ( ed_Ilefts - ed_Irights )
# print(ed_Ilefts)
# print(ed_Irights)
# print(ed_Imean)

print(ed_Imean_old)
print(ed_Imean)

dt = 0.014
ts = linspace(4.2, 6.3, num=16)
order = 7
k = 5
chi = 160
# our data
U = 0.

Ilefts_0_old, Irights_0_old, bd_0_old = read_tempo_data(Vs, U)
Imean_0_old = 0.5 * ( Ilefts_0_old - Irights_0_old )

Ilefts_0, Irights_0, bd_0 = read_tempo_data2(Vs, U, dt=dt, order=order, k=k, chi=chi)
Imean_0 = [0.5 * ( a - b ) for a, b in zip(Ilefts_0, Irights_0)]

# print(Imean_0[4])

# print(Ilefts_0[4])
# print(Irights_0[4])

U = 2.
Ilefts_2_old, Irights_2_old, bd_2_old = read_tempo_data(Vs, U)
Imean_2_old = 0.5 * ( Ilefts_2_old - Irights_2_old )


Ilefts_2, Irights_2, bd_2 = read_tempo_data2(Vs, U, dt=dt, order=order, k=k, chi=chi)
Imean_2 = [0.5 * ( a - b ) for a, b in zip(Ilefts_2, Irights_2)]

# print([item[0] for item in Imean_2])
# print([item[-1] for item in Imean_2])

U = 4.
Ilefts_4_old, Irights_4_old, bd_4_old = read_tempo_data(Vs, U)
Imean_4_old = 0.5 * ( Ilefts_4_old - Irights_4_old )

Ilefts_4, Irights_4, bd_4 = read_tempo_data2(Vs, U, dt=dt, order=order, k=k, chi=chi)
Imean_4 = [0.5 * ( a - b ) for a, b in zip(Ilefts_4, Irights_4)]

# print(Imean_4_old[0])
# print(Imean_4[0])

# print([item[0] for item in Imean_4])
# print([item[-1] for item in Imean_4])

# U = 6.
# Ilefts_6, Irights_6, bd_6 = read_tempo_data2(Vs, U, dt=dt, order=order, k=k, chi=chi)
# Imean_6 = [0.5 * ( a - b ) for a, b in zip(Ilefts_6, Irights_6)]

# U = 8.
# Ilefts_8, Irights_8, bd_8 = read_tempo_data(Vs2, U)
# Imean_8 = 0.5 * ( Ilefts_8 - Irights_8 )


# print(Ilefts_8)
# print(Irights_8)
# print(asarray(Ilefts_8) + asarray(Irights_8))

Imeans_allU_0 = [Imean_0[0], Imean_2[0], Imean_4[0]]
Imeans_allU_1 = [Imean_0[1], Imean_2[1], Imean_4[1]]
Imeans_allU_2 = [Imean_0[2], Imean_2[2], Imean_4[2]]
Imeans_allU_3 = [Imean_0[3], Imean_2[3], Imean_4[3]]

Imeans_allU_old_0 = [Imean_0_old[0], Imean_2_old[0], Imean_4_old[0]]
Imeans_allU_old_1 = [Imean_0_old[1], Imean_2_old[1], Imean_4_old[1]]
Imeans_allU_old_2 = [Imean_0_old[2], Imean_2_old[2], Imean_4_old[2]]
Imeans_allU_old_3 = [Imean_0_old[3], Imean_2_old[3], Imean_4_old[3]]


fontsize = 16
labelsize = 14
linewidth = 1.
markersize = 6

fig, ax = plt.subplots(2, 2, figsize=(8, 6))


ed_ts = [ts[0], ts[-1]]
ts_old = [ts[0]]

Us = [0,2,4]
# U = 0

ed_I = asarray([ed_Imean_old[0], ed_Imean[0]])
ax[0,0].plot(ed_ts, ed_I , linewidth=1.5, color='r', marker=markers[0], markersize=markersize, markerfacecolor='none', ls = 'none')
for i, U in enumerate(Us):
	# ed_I = asarray([ed_Imean_old[i], ed_Imean[i]])
	# ax[0,0].plot(ed_ts, ed_I , linewidth=1.5, color='r', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = 'none')

	ax[0,0].plot(ts_old, [Imeans_allU_old_0[i]] , linewidth=1.5, color='c', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = 'none')
	ax[0,0].plot(ts, Imeans_allU_0[i] , linewidth=1.5, color='g', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = '--', label=r'$U/\Gamma=%s$'%(U))

# ax.plot(ts, Imean_0[4], linewidth=1.5, color=colors[0],  markersize=markersize, markerfacecolor='none', ls = '--')

# ax.text(2.,1.7,r'$U/\Gamma = 0$',rotation=30, fontsize=10)
# ax.text(2.,1.48,r'$U/\Gamma = 2$',rotation=30, fontsize=10)
# ax.text(2.,1.08,r'$U/\Gamma = 4$',rotation=25, fontsize=10)
# ax.text(2.,0.75,r'$U/\Gamma = 6$',rotation=20, fontsize=10)
# ax.text(2.,0.54,r'$U/\Gamma = 8$',rotation=15, fontsize=10)


ax[0,0].set_ylabel(r'$2\pi \mathcal{J}/V$', fontsize=fontsize)
ax[0,0].set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax[0,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,0].tick_params(axis='y', which='both')
ax[0,0].locator_params(axis='both', nbins=6)
# ax.set_ylim(0, 2.)
# ax.set_xlim(0, 3)
ax[0,0].set_title(r'$V/\Gamma=%s$'%(Vs[0]), fontsize=fontsize)
ax[0,0].legend(fontsize=10)



ed_I = asarray([ed_Imean_old[1], ed_Imean[1]])
# print(ed_I)
ax[0,1].plot(ed_ts, ed_I , linewidth=1.5, color='r', marker=markers[0], markersize=markersize, markerfacecolor='none', ls = 'none')
for i, U in enumerate(Us):

	ax[0,1].plot(ts_old, [Imeans_allU_old_1[i]] , linewidth=1.5, color='c', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = 'none')
	ax[0,1].plot(ts, Imeans_allU_1[i] , linewidth=1.5, color='g', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = '--', label=r'$U/\Gamma=%s$'%(U))


ax[0,1].set_ylabel(r'$2\pi \mathcal{J}/V$', fontsize=fontsize)
ax[0,1].set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax[0,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0,1].tick_params(axis='y', which='both')
ax[0,1].locator_params(axis='both', nbins=6)
ax[0,1].set_title(r'$V/\Gamma=%s$'%(Vs[1]), fontsize=fontsize)

ed_I = asarray([ed_Imean_old[2], ed_Imean[2]])
# print(ed_I)
ax[1,0].plot(ed_ts, ed_I , linewidth=1.5, color='r', marker=markers[0], markersize=markersize, markerfacecolor='none', ls = 'none')
for i, U in enumerate(Us):

	ax[1,0].plot(ts_old, [Imeans_allU_old_2[i]] , linewidth=1.5, color='c', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = 'none')
	ax[1,0].plot(ts, Imeans_allU_2[i] , linewidth=1.5, color='g', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = '--', label=r'$U/\Gamma=%s$'%(U))


ax[1,0].set_ylabel(r'$2\pi \mathcal{J}/V$', fontsize=fontsize)
ax[1,0].set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax[1,0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,0].tick_params(axis='y', which='both')
ax[1,0].locator_params(axis='both', nbins=6)
ax[1,0].set_title(r'$V/\Gamma=%s$'%(Vs[2]), fontsize=fontsize)

ed_I = asarray([ed_Imean_old[3], ed_Imean[3]])
# print(ed_I)
ax[1,1].plot(ed_ts, ed_I , linewidth=1.5, color='r', marker=markers[0], markersize=markersize, markerfacecolor='none', ls = 'none')
for i, U in enumerate(Us):

	ax[1,1].plot(ts_old, [Imeans_allU_old_3[i]] , linewidth=1.5, color='c', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = 'none')
	ax[1,1].plot(ts, Imeans_allU_3[i] , linewidth=1.5, color='g', marker=markers[i], markersize=markersize, markerfacecolor='none', ls = '--', label=r'$U/\Gamma=%s$'%(U))


ax[1,1].set_ylabel(r'$2\pi \mathcal{J}/V$', fontsize=fontsize)
ax[1,1].set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax[1,1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1,1].tick_params(axis='y', which='both')
ax[1,1].locator_params(axis='both', nbins=6)
ax[1,1].set_title(r'$V/\Gamma=%s$'%(Vs[3]), fontsize=fontsize)

plt.tight_layout(pad=0.5)

# plt.savefig('fig3.pdf', bbox_inches='tight')

plt.show()

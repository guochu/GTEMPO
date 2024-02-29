from matplotlib import pyplot as plt
from matplotlib import rc

import pandas
import json
from numpy import asarray

rc('text', usetex=True)

def read_tempo_data(t, t0, U, dt, order, chi):
	mu = U / 2
	filename = 'result/anderson_tempo1_beta40.0_t%s_%s_U%s_mu%s_dt%s_order%s_chi%s.json'%(t, t0, U, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return asarray(data['ts']) * 0.1, asarray(data['nn']), asarray(data['nup']), asarray(data['ndown']), asarray(data['vacuum']), asarray(data['bd'])

def read_tempo_data_2(t, t0, U, dt, order, chi):
	mu = U / 2
	filename = 'result/anderson_tempo1_beta40.0_t%s_%s_U%s_mu%s_dt%s_order%s_chi%s_2.json'%(t, t0, U, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	return asarray(data['ts']) * 0.1, asarray(data['nn']), asarray(data['nup']), asarray(data['ndown']), asarray(data['vacuum']), asarray(data['bd'])



fontsize = 16
labelsize = 14
linewidth = 2.
markersize = 6

colors = ['b', 'g', 'c', 'y', 'r']

fig, ax = plt.subplots(1, 2, figsize=(8, 4.))

t = 25.
dt = 0.05
t0 = 5.
U = 1.
order = 10
chi = 60

tempo_data = read_tempo_data(t=t, t0=t0, U=U, dt=dt, order=order, chi=chi)


ax[0].plot(tempo_data[0], tempo_data[1], linewidth=linewidth, ls='--', color=colors[0], label=r'$\vert \uparrow\downarrow\rangle$')
ax[0].plot(tempo_data[0], tempo_data[2], linewidth=linewidth, ls='--', color=colors[1], label=r'$\vert \uparrow\rangle$')
ax[0].plot(tempo_data[0], tempo_data[3], linewidth=linewidth, ls='--', color=colors[2], label=r'$\vert \downarrow\rangle$')
ax[0].plot(tempo_data[0], tempo_data[4], linewidth=linewidth, ls='--', color=colors[3], label=r'$\vert 0\rangle$')


ax[0].set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax[0].set_ylabel(r'Populations', fontsize=fontsize)
ax[0].tick_params(axis='both', which='major', labelsize=labelsize)
ax[0].tick_params(axis='y', which='both')
ax[0].locator_params(axis='both', nbins=6)
ax[0].set_ylim(0, 1)
ax[0].set_xlim(0, t * 0.1)
ax[0].annotate(r'(a)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)
ax[0].legend(fontsize=12)


tempo_data = read_tempo_data_2(t=t, t0=t0, U=U, dt=dt, order=order, chi=chi)

ax[1].plot(tempo_data[0], tempo_data[1], linewidth=linewidth, ls='--', color=colors[0], label=r'$\vert \uparrow\downarrow\rangle$')
ax[1].plot(tempo_data[0], tempo_data[2], linewidth=linewidth, ls='--', color=colors[1], label=r'$\vert \uparrow\rangle$')
ax[1].plot(tempo_data[0], tempo_data[3], linewidth=linewidth, ls='--', color=colors[2], label=r'$\vert \downarrow\rangle$')
ax[1].plot(tempo_data[0], tempo_data[4], linewidth=linewidth, ls='--', color=colors[3], label=r'$\vert 0\rangle$')


ax[1].set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax[1].set_ylabel(r'Populations', fontsize=fontsize)
ax[1].tick_params(axis='both', which='major', labelsize=labelsize)
ax[1].tick_params(axis='y', which='both')
ax[1].locator_params(axis='both', nbins=6)
# ax[1].legend(fontsize=12)
ax[1].set_ylim(0, 1)
ax[1].set_xlim(0, t * 0.1)
ax[1].annotate(r'(b)', xy=(0.1, 0.85),xycoords='axes fraction', fontsize=fontsize)




plt.tight_layout(pad=0.5)


# plt.savefig('fig3.pdf', bbox_inches='tight')

plt.show()
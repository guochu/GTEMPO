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

def read_mixed_tempo(beta, t, U, dt=0.05, dtau=0.1, order=10, chi=60):
	# mpath = '/Users/guochu/Documents/Missile/GTEMPO/examples/thirdpaper/oneorbital'
	mu = U/2
	# t2 = t/2
	filename = 'result/mixed_beta%s_dtau%s_t%s_dt%s_U%s_mu%s_order%s_chi%s.json'%(beta, dtau, t, dt, U, mu, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	nup = data['nup']
	ndown = data['ndown']
	nn = data['nn']
	vacuum = data['vacuum']
	ts = asarray(data['ts'])[:-1]
	return ts-ts[0], -asarray(nup), -asarray(ndown), asarray(nn), vacuum


fontsize = 20
labelsize = 16
linewidth = 2.5
markersize = 4


colors = ['b', 'c', 'orange', 'g', 'y', 'r', 'k']
markers = ['o', '^', '+']
color0 = 'g'
color1 = 'r'

fig, ax = plt.subplots(1, 1, figsize=(6,5))


beta = 40.
dt = 0.05

U = 0.5

t0 = 20.

chi_r = 100

dtau = 0.1


order = 10

k = 8
	
chi_m = 80
mixed_ts, nup, ndown, nn, vacuum = read_mixed_tempo(beta, t0, U, dt, dtau=dtau, chi=chi_m)
ns = nn * 2 + nup + ndown

ax.plot(mixed_ts * 0.1, ns-1, linewidth=linewidth, ls='--', color=colors[0], label=r'$\chi=%s$'%(chi_m))

chi_m = 120
mixed_ts, nup, ndown, nn, vacuum = read_mixed_tempo(beta, t0, U, dt, dtau=dtau, chi=chi_m)
ns = nn * 2 + nup + ndown

ax.plot(mixed_ts * 0.1, ns-1, linewidth=linewidth, ls='--', color=colors[1], label=r'$\chi=%s$'%(chi_m))


# ax.legend(fontsize=10)


ax.set_xlabel(r'$\Gamma t$', fontsize=fontsize)
ax.set_ylabel(r'$\bar{n}-1$', fontsize=fontsize)
ax.tick_params(axis='both', which='major', labelsize=labelsize)
ax.locator_params(axis='both', nbins=6)
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
ax.legend(fontsize=12)




plt.tight_layout(pad=0.5)

plt.savefig('oneorbital2.pdf', bbox_inches='tight')

plt.show()


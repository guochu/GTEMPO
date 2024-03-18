import json
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.colors import LogNorm
from numpy import linspace, asarray, loadtxt, linspace, sqrt
from numpy.linalg import norm
import math


def parse_complex_array(data):
	re = [item['re'] for item in data]
	im = [item['im'] for item in data]
	return asarray(re) + 1j * asarray(im)


def read_real_tempo_2(beta, U, dt, order=10, chi=60):
	t0 = 60.
	t = 100.
	mu = U/2
	# t2 = t/2
	filename = 'result/anderson_tempo1_beta%s_t%s_%s_U%s_mu%s_dt%s_order%s_chi%s.json'%(beta, t, t0, U, mu, dt, order, chi)
	with open(filename, 'r') as f:
		data = f.read()
		data = json.loads(data)
	gt = data['gt'][:401]
	lt = data['lt'][:401]

	gf_ts = data['gf_ts'][:401]

	new_data = {'gt':gt, 'lt':lt, 'gf_ts':gf_ts, 'ns':data['ns'], 'ts':data['ts'], 'bd':data['bd']}

	filename2 = 'result/anderson_tempo1_beta%s_t%s_%s_U%s_mu%s_dt%s_order%s_chi%s.json'%(beta, 80., t0, U, mu, dt, order, chi)

	with open(filename2, 'w') as f:
		json.dump(new_data, f)




beta = 40.
U = 0.1
dt = 0.05

for U in [0.1, 0.5, 1.]:
	x = read_real_tempo_2(beta, U, dt)
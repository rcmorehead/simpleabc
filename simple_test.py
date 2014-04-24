
import matplotlib
matplotlib.use('Agg')
import ABC 
import simplest_model  
#from astropy.io import ascii
import numpy as np
import pickle 
from pylab import * 
import time 


def main():
	#stars = ascii.read('stars.csv')
	stars = pickle.load(file('stars.pkl'))
	#obs = ascii.read('test_cat.dat')
	obs = pickle.load(file('data.pkl'))

	out = file('benchmark.dat', 'w')

	eps = 0.01
	min_part = 100

	model = simplest_model.MyModel(stars, obs)

	n_procs = [1,2,3,4,5,6,7,8]

	start = time.time()
	post, reject, accept, tot = ABC.basic_abc(model, obs, min_particles=min_part, 
												target_epsilon=eps, parallel=False)
	end = time.time()
	print 'Serial took {}s'.format(end - start)
	print >> out, '0,{},{}'.format(end - start,(end - start)/float(tot))

	re_bi = [x[0] for x in reject]
	re_sc = [x[1] for x in reject]
	ac_bi = [x[0] for x in post]
	ac_sc = [x[1] for x in post]

	figure()
	ax=subplot(111)
	plot(re_bi, re_sc, 'o', c='0.5',mec='0.5')
	axvline(5, ls='--',color='orange')
	axhline(3, ls='--',color='orange')
	plot(ac_bi,ac_sc,'o', c='k')
	title('Serial total:{} accepted:{} epsilon={}'.format(tot, accept, eps))
	xlabel('n')
	ylabel('sigma')
	savefig('Serial')

	for N in n_procs:

		start = time.time()
		post, reject, accept, tot = ABC.basic_abc(model, obs, min_particles=min_part, 
													target_epsilon=eps, parallel=True, n_procs=N)
		end = time.time()
		print 'Parallel took {}s'.format(end - start)
		print >> out, '{},{},{}'.format(N,end - start,(end - start)/float(tot))

		re_bi = [x[0] for x in reject]
		re_sc = [x[1] for x in reject]
		ac_bi = [x[0] for x in post]
		ac_sc = [x[1] for x in post]

		figure()
		ax=subplot(111)
		plot(re_bi, re_sc, 'o', c='0.5',mec='0.5')
		axvline(5, ls='--',color='orange')
		axhline(3, ls='--',color='orange')
		plot(ac_bi,ac_sc,'o', c='k')
		xlabel('n')
		ylabel('sigma')
		title('Parallel total:{} accepted:{} epsilon={}'.format(tot, accept, eps))
		savefig('Parallel_{}_procs'.format(N))


	out.close()


if __name__ == '__main__':
	main()


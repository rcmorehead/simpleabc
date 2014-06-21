
import matplotlib
matplotlib.use('Agg')
from simple_abc import *
import simple_model
#from astropy.io import ascii
import numpy as np
import pickle 
from pylab import * 
import time 
from scipy import stats
from triangle import corner


def main():
    np.random.seed(917)

    steps = 5
    eps = 0.25
    min_part = 10

    #stars = pickle.load(file('stars.pkl'))
    stars = pickle.load(file('stars_trimmed.pkl'))
    #obs = pickle.load(file('data.pkl'))

    model = simple_model.MyModel(stars)
    model.set_prior([stats.uniform(0.5, 1.0),
                    stats.uniform(0, 1.0)])

    theta = (0.513265306122, 0.1)

    obs = model.generate_data(theta)
    model.set_data(obs)





    n_procs = [1, 2, 3, 4, 5, 6, 7, 8]

    start = time.time()
    OT = pmc_abc(model, obs, epsilon_0=eps, min_particles=min_part, steps=steps,
    						 target_epsilon=eps, parallel=False, plot=True)
    end = time.time()
    print 'Serial took {}s'.format(end - start)
    out_pickle = file('simptest.pkl', 'w')
    pickle.dump(OT, out_pickle)
    out_pickle.close()
   # for step in xrange(steps):
   #     figure()
   #     corner(OT[step]['theta accepted'], truths=(theta[0], theta[1]))
   #     #show()
   #     savefig('PLOTS/simple{}_posterior.png'.format(step))
    #     THETAR = asarray(OT[step]['theta rejected']).T
    #     THETAA = asarray(OT[step]['theta accepted']).T
    #     figure()
    #     plot(THETAR[0],
    #              THETAR[1], 'o', c='k',
    #              alpha=.1)
    #
    #     plot(THETAA[0],
    #              THETAA[1] , 'o', mfc='None', mec='m', mew=1)
    #
    #     axvline(mu, linestyle='--',lw=3)
    #     xlabel(r'$\mu$', fontsize=24)
    #     ylabel(r'$\sigma$', fontsize=24)
    #     axhline(sigma, linestyle='--',lw=3)
    #     title("Step {}".format(step), fontsize=24)
    #     xlim(-2,2)
    #     ylim(0,4)
    #     savefig('PLOTS/{}.png'.format(step))
    #
    #     x1 = linspace(THETAA[0].min(), THETAA[0].max())
    #     x2 = linspace(THETAA[1].min(), THETAA[1].max())
    #
    #     figure()
    #     subplot(121)
    #     hist(THETAA[0], bins=sqrt(THETAA[0].size), alpha=.2, fc='gray',
    #          normed=True)
    #     gkde = stats.gaussian_kde(THETAA[0])
    #     plot(x1,gkde(x1), lw=2)
    #     axvline(mu, lw=2)
    #     xlabel(r'$\mu$', fontsize=24)
    #     subplot(122)
    #     hist(THETAA[1], bins=sqrt(THETAA[1].size), alpha=.2, fc='gray',
    #          normed=True)
    #     gkde = stats.gaussian_kde(THETAA[1])
    #     plot(x2,gkde(x2), lw=2)
    #     axvline(sigma, lw=2)
    #     xlabel(r'$\sigma$', fontsize=24)
    #     suptitle("Step {}".format(step), fontsize=24)
    #     savefig('PLOTS/hist_{}.png'.format(step))
    # figure()
    # plot(range(steps), OT[:]['epsilon'], 'o-')
    # title(r'$\epsilon$', fontsize=24)
    # ylabel(r'$\epsilon$', fontsize=24)
    # xlabel('step', fontsize=24)
    # savefig('PLOTS/epsilon.png')


if __name__ == '__main__':
	main()


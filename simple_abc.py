"""
Module for Approximate Bayesian Computation

"""

import multiprocessing
import numpy as np
from scipy import stats

def basic_abc(model, data, target_epsilon=0.1, min_particles=10,
              parallel=False, n_procs='all'):
    """
    Preform Approximate Bayesian Computation on a data set given a forward
    model.

    """


    posterior, rejected = [], []
    epsilon = target_epsilon
    trial_count, accepted_count = 0, 0

    data_summary_stats = model.summary_stats(data)

    if parallel:
        attempts = 2*min_particles
        if n_procs == 'all':
            n_procs = multiprocessing.cpu_count()
        while accepted_count < min_particles:
            thetas = [model.draw_theta() for x in
                               xrange(attempts)]

            #Start a pool of workers



            pool = multiprocessing.Pool(n_procs)
            ds = pool.map(model, thetas)

            #Shut down pool
            pool.close()
            pool.join()

            for j, d in enumerate(ds):
                if d < epsilon:
                    posterior.append(thetas[j])
                    accepted_count += 1
                    trial_count += 1
                else:
                    rejected.append(thetas[j])
                    trial_count += 1

            attempts = int(float(trial_count)/float(accepted_count + 1) *
                        (min_particles - accepted_count))

        return posterior, rejected, accepted_count, trial_count
    else:
        while accepted_count <= min_particles:

            trial_count += 1


            theta = model.draw_theta()
            synthetic_data = model.generate_data(theta)

            synthetic_summary_stats = model.summary_stats(synthetic_data)
            distance = model.distance_function(data_summary_stats,
                                               synthetic_summary_stats)

            if distance < epsilon:
                accepted_count += 1
                posterior.append(theta)
            else:
                rejected.append(theta)

        return posterior, rejected, accepted_count, trial_count


def pmc_abc(model, data, target_epsilon=0.1, epsilon_0=0.25, min_particles=10,
              parallel=False, n_procs='all'):
    """
    Preform Approximate Bayesian Computation on a data set given a forward
    model using pmc.

    """
    #Fist ABC calculation
    posterior, rejected, accepted_count, trial_count = basic_abc(model,
                                                                 data,
                                                                 epsilon_0,
                                                                 min_particles,
                                                                 parallel,
                                                                 n_procs)
    #TODO make theta-like vectors constant size when generated
    theta = np.asarray(posterior)
    theta = theta.T
    tau2 = np.zeros((1, theta.shape[0]))

    weights = np.ones((theta.shape[0], theta[1].size))


    theta_star = np.zeros_like(theta)
    theta_i = np.zeros_like(theta)
    for j in xrange(theta.shape[0]):
        tau2[0][j] = 2*np.var(theta[j])

        weights[j] = weights[j]*1/float(theta[j].size)

        theta_star[j] = np.random.choice(theta[j],
                                         size=theta[j].size, replace=True,
                                         p=weights[j])
        theta_i[j] = stats.norm.rvs(loc=theta_star[j], scale=np.sqrt(tau2[0][j]))


    #print theta_i
    new_weights = calc_weights(theta_i, theta, tau2, weights,
                               prior=stats.uniform(-2, 4))
                               #TODO add prior specifcation tp model
    print new_weights
    return posterior, rejected, accepted_count, trial_count


def calc_weights(theta_i, theta, tau2, weights,prior="None"):

    """
    Calculates importance weights
    """
    weights_new = np.zeros_like(theta_i)

    for i in xrange(theta_i.shape[0]):
        for j in xrange(theta_i[i].size):
            print theta_i[i][j], prior.pdf(theta_i[i][j]),np.sum(weights[i]*stats.norm.pdf((theta_i[i][j] -theta[i])/np.sqrt(tau2[0][i])))
            weights_new[i][j] = (prior.pdf(theta_i[i][j]) /
                              np.sum(weights[i]*stats.norm.pdf((theta_i[i][j] - theta)/np.sqrt(tau2[0][i]))))

    return weights_new
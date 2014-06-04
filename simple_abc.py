"""
Module for Approximate Bayesian Computation

"""

import multiprocessing
import numpy as np

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


    tau2 = 2*np.var(posterior)

    weights = np.ones(posterior.size).fill(1/posterior.size)

    theta_star = np.random.choice(posterior, size=min_particles*10,
                                  replace=True, p=weights)


    print weights,theta_star
    return posterior, rejected, accepted_count, trial_count
"""
Module for Approximate Bayesian Computation

"""

import multiprocessing


def basic_abc(model, data, target_epsilon=0.1, min_particles=10,
              parallel=False):
    """
    Preform Approximate Bayesian Computation on a data set given a forward
    model.

    """


    posterior, rejected = [], []
    epsilon = target_epsilon
    trial_count, accepted_count = 0, 0

    data_summary_stats = model.summary_stats(data)

    if parallel:
        number_of_particles = min_particles * 2
        while accepted_count < min_particles:
            thetas = [model.draw_theta() for x in xrange(number_of_particles)]
            print thetas
            #Start a pool of workers
            N_procs = multiprocessing.cpu_count()
            print N_procs
            N_procs = 1
            pool = multiprocessing.Pool(N_procs)

            sum_stats = pool.map(model, thetas)

            #Shut down pool
            pool.close()
            pool.join()

            print sum_stats
            accepted_count = min_particles

        return 'This ', 'Just ', 'Worked ', 'Yo'
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
"""
Module for Approximate Bayesian Computation

Usage
-----

    ABC.py <file containing stars> <model.py> <observed>


"""


def abc(model, data, target_epsilon=0.1, min_particles=10, parallel=False):
    """
    Preform Approximate Bayesian Computation on a data set given a forward
    model.

    """
    if parallel:
        raise NotImplementedError('Parallelism coming soon. I promise. ')

    posterior, rejected = [], []
    epsilon = target_epsilon
    trial_count, accepted_count = 0, 0

    data_summary_stats = model.summary_stats(data)

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
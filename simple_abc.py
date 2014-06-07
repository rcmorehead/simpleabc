"""
Module for Approximate Bayesian Computation

"""
from abc import ABCMeta, abstractmethod
import multiprocessing
import numpy as np
from scipy import stats


class Model(object):
    """
    Base class for constructing models for approximate bayesian computing
    and various uses limited only by the user's imagination.

    WARNING!! Not meant for direct use! You must implement your own model as a
    subclass and override the all following methods:

    * Model.draw_theta
    * Model.generate_data
    * Model.summary_stats
    * Model.distance_function

    """

    __metaclass__ = ABCMeta

    def __call__(self, theta):
        return self.generate_data_and_reduce(theta)

    def set_data(self, data):
        self.data = data
        self.data_sum_stats = self.summary_stats(self.data)

    #TODO think about a beter way to handle prior functions
    def set_prior(self, prior):
        self.prior = prior

    def generate_data_and_reduce(self, theta):
        """
        A combined method for generating data, calculating summary statistics
        and evaluating the distance function all at once.
        """
        synth = self.generate_data(theta)
        sum_stats = self.summary_stats(synth)
        d = self.distance_function(sum_stats, self.data_sum_stats)

        return d


    @abstractmethod
    def draw_theta(self):
        """
        Sub-classable method for drawing from a prior distribution.

        This method should return an array-like iterable that is a vector of
        proposed model parameters from your prior distribution.
        """

    @abstractmethod
    def generate_data(self, theta):
        """
        Sub-classable method for generating synthetic data sets from forward
        model.

        This method should return an array/matrix/table of simulated data
        taking vector theta as an argument.
        """

    @abstractmethod
    def summary_stats(self, data):
        """
        Sub-classable method for computing summary statistics.

        This method should return an array-like iterable of summary statistics
        taking an array/matrix/table as an argument.
        """

    @abstractmethod
    def distance_function(self, summary_stats, summary_stats_synth):
        """
        Sub-classable method for computing a distance function.

        This method should return a distance D of for comparing to the
        acceptance tolerance (epsilon) taking two array-like iterables of
        summary statistics as an argument (nominally the observed summary
        statistics and .
        """
################################################################################
#########################    ABC Algorithms   ##################################
################################################################################

def basic_abc(model, data, epsilon=0.1, min_particles=10,
              parallel=False, n_procs='all', pmc_mode=False,
              weights='None', theta_prev='None', tau_squared='None'):
    """
    Preform Approximate Bayesian Computation on a data set given a forward
    model.

    """


    posterior, rejected, distances = [], [], []
    trial_count, accepted_count = 0, 0

    data_summary_stats = model.summary_stats(data)
    #TODO Implement pmc option in parallel mode
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

        return (posterior, rejected, distances,
                accepted_count, trial_count,
                epsilon)
    else:
        while accepted_count <= min_particles:

            trial_count += 1

            if pmc_mode:
                theta_star = list(theta_prev.shape[0])
                theta = list(theta_prev.shape[0])

                for j in xrange(theta_prev.shape[0]):
                    theta_star[j] = np.random.choice(theta_prev[j],
                                         replace=True,
                                         p=weights[j])
                    theta[j] = stats.norm.rvs(loc=theta_star[j],
                                    scale=np.sqrt(tau_squared[0][j]))

            else:
                theta = model.draw_theta()

            synthetic_data = model.generate_data(theta)

            synthetic_summary_stats = model.summary_stats(synthetic_data)
            distance = model.distance_function(data_summary_stats,
                                               synthetic_summary_stats)

            if distance < epsilon:
                accepted_count += 1
                posterior.append(theta)
                distances.append(distance)
            else:
                rejected.append(theta)

        return (posterior, rejected, distances,
                accepted_count, trial_count,
                epsilon)


def pmc_abc(model, data, target_epsilon=0.1, epsilon_0=0.25, min_particles=1000,
              steps=10, parallel=False, n_procs='all'):
    """
    Preform Approximate Bayesian Computation on a data set given a forward
    model using pmc.

    """
    output_record = np.empty(steps, dtype=[('theta accepted', object),
                                           ('theta rejected', object),
                                           ('D accepted', object),
                                           ('n accepted', float),
                                           ('n total', float),
                                           ('epsilion', float),
                                           ] )

    epsilon = epsilon_0

    for step in xrange(steps):
        if step == 0:
    #Fist ABC calculation
            output_record[step] = basic_abc(model, data, epsilon,
                                            min_particles, parallel,
                                            n_procs, pmc_mode=False)

            theta_prev = np.asarray(output_record[step]['theta accepted']).T
            tau_squared = np.zeros((1, theta_prev.shape[0]))
            weights = np.ones((theta_prev.shape[0], theta_prev[1].size))

            for j in xrange(theta_prev.shape[0]):
                tau_squared[0][j] = 2*np.var(theta_prev[j])
                weights[j] = weights[j]*1/float(theta_prev[j].size)

            epsilon = stats.scoreatpercentile(output_record[step]['D accepted'],
                                              per=75)
            print tau_squared
            print weights
            print epsilon

        else:
            pass


        # theta = np.asarray(posterior)
        # theta = theta.T
        # tau_squared = np.zeros((1, theta.shape[0]))
        #
        # weights = np.ones((theta.shape[0], theta[1].size))
        #
        # theta_i, rej, acc, trial = basic_abc(model, data, target_epsilon=0.1,
        #                                  min_particles=10, parallel=False,
        #                                  n_procs='all', pmc_mode=False,
        #                                  weights=weights, theta_prev=theta,
        #                                  tau_squared=tau_squared)
        #
        #
        # theta_i = np.asarray(theta_i).T
        #
        # new_weights = calc_weights(theta_i, theta, tau_squared, weights,
        #                        prior=stats.uniform(-2, 4))

        #mean.prev <- sum(mu[t-1,]*weights[t-1,])
		#var.prev <- sum((mu[t-1,] - mean.prev)^2*weights[t-1,])
        #
        #                       #TODO add prior specifcation tp model

    print output_record

    return output_record


def calc_weights(theta_i, theta, tau_squared, weights,prior="None"):

    """
    Calculates importance weights
    """
    weights_new = np.zeros_like(theta_i)

    for i in xrange(theta_i.shape[0]):
        for j in xrange(theta_i[i].size):
            weights_new[i][j] = (prior.pdf(theta_i[i][j]) /
                                np.sum(weights[i]*stats.norm.pdf(
                                  (theta_i[i][j] - theta),
                                  np.sqrt(tau_squared[0][i]))))


    return weights_new
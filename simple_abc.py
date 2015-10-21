"""
Module for Approximate Bayesian Computation

"""
from abc import ABCMeta, abstractmethod
import multiprocessing
import numpy as np
from scipy import stats
from numpy.lib.recfunctions import stack_arrays
from numpy.testing import assert_almost_equal

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

def basic_abc(model, data, epsilon=1, min_samples=10,
              parallel=False, n_procs='all', pmc_mode=False,
              weights='None', theta_prev='None', tau_squared='None'):
    """
    Perform Approximate Bayesian Computation (ABC) on a data set given a
    forward model.

    ABC is a likelihood-free method of Bayesian inference that uses simulation
    to approximate the true posterior distribution of a parameter. It is
    appropriate to use in situations where:

    The likelihood function is unknown or is too computationally
    expensive to compute.

    There exists a good forward model that can produce data sets
    like the one of interest.

    It is not a replacement for other methods when a likelihood
    function is available!

    Parameters
    ----------
    model : object
        A model that is a subclass of simpleabc.Model
    data  : object, array_like
        The "observed" data set for inference.
    epsilon : float, optional
        The tolerance to accept parameter draws, default is 1.
    min_samples : int, optional
        Minimum number of posterior samples.
    parallel : bool, optional
        Run in parallel mode. Default is a single thread.
    n_procs : int, str, optional
        Number of subprocesses in parallel mode. Default is 'all' one for each
        available core.
    pmc_mode : bool, optional
        Population Monte Carlo mode on or off. Default is False. This is not
        meant to be called by the user, but is set by simple_abc.pmc_abc.
    weights : object, array_like, str, optional
        Importance sampling weights from previous PMC step. Used  by
        simple_abc.pmc_abc only.
    theta_prev : object, array_like, str, optional
        Posterior draws from previous PMC step.  Used by simple_abc.pmc_abc
        only.
    tau_squared : object, array_like, str, optional
        Previous Gaussian kernel variances. for importance sampling. Used by
        simple_abc.pmc_abc only.

    Returns
    -------
    posterior : numpy array
        Array of posterior samples.
    distances : object
        Array of accepted distances.
    accepted_count : float
        Number of  posterior samples.
    trial_count : float
        Number of total samples attempted.
    epsilon : float
        Distance tolerance used.
    weights : numpy array
        Importance sampling weights. Returns an array of 1s where
        size = posterior.size when not in pmc mode.
    tau_squared : numpy array
        Gaussian kernel variances. Returns an array of 0s where
        size = posterior.size when not in pmc mode.
    eff_sample : numpy array
        Effective sample size. Returns an array of 1s where
        size = posterior.size when not in pmc mode.

    Examples
    --------
    Forth coming.

    """


    posterior, rejected, distances = [], [], []
    trial_count, accepted_count = 0, 0


    data_summary_stats = model.summary_stats(data)
    #TODO Implement pmc option in parallel mode
    if parallel:
        attempts = 2*min_samples 
        if n_procs == 'all':
            n_procs = multiprocessing.cpu_count()
        while accepted_count < min_samples :
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
                    #rejected.append(thetas[j])
                    trial_count += 1

            attempts = int(float(trial_count)/float(accepted_count + 1) *
                        (min_samples - accepted_count))

        return (posterior, distances,
                accepted_count, trial_count,
                epsilon)
    else:
        while accepted_count <= min_samples:

            trial_count += 1

            if pmc_mode:
                #theta_star = []
                #theta = []

                #for j in xrange(theta_prev.shape[0]):
                #    theta_star.append(np.random.choice(theta_prev[j],
                #                         replace=True,
                #                         p=weights[j]))
                #    #print "t*,tu2: ",theta_star[j], np.sqrt(tau_squared[0][j])
                #    theta.append(stats.norm.rvs(loc=theta_star[j],
                #                    scale=np.sqrt(tau_squared[0][j])))
                #print theta_prev
                theta_star = theta_prev[:, np.random.choice(
                                        xrange(0, theta_prev.shape[1]),
                                        replace=True, p=weights/weights.sum())]

                theta = stats.multivariate_normal.rvs(theta_star, tau_squared)
                if np.isscalar(theta) == True:
                    theta = [theta]


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
                pass
                #rejected.append(theta)

        posterior = np.asarray(posterior).T

        if len(posterior.shape) > 1:
            n = posterior.shape[1]
        else:
            n = posterior.shape[0]


        weights = np.ones(n)
        tau_squared = np.zeros((posterior.shape[0], posterior.shape[0]))
        eff_sample = n

        return (posterior, distances,
                accepted_count, trial_count,
                epsilon, weights, tau_squared, eff_sample)


def pmc_abc(model, data, epsilon_0=1, min_samples=10,
              steps=10, resume=None, parallel=False, n_procs='all'):
    """
    Perform a sequence of ABC posterior approximations using the sequential
    population Monte Carlo algorithm.


    Parameters
    ----------
    model : object
        A model that is a subclass of simpleabc.Model
    data  : object, array_like
        The "observed" data set for inference.
    epsilon_0 : float, optional
        The initial tolerance to accept parameter draws, default is 1.
    min_samples : int, optional
        Minimum number of posterior samples.
    steps : int
        The number of pmc steps to attempt
    resume : numpy record array, optional
        A record array of a previous pmc sequence to continue the sequence on.
    parallel : bool, optional
        Run in parallel mode. Default is a single thread.
    n_procs : int, str, optional
        Number of subprocesses in parallel mode. Default is 'all' one for each
        available core.

    Returns
    -------

    output_record : numpy record array
        A record array containing all ABC output for each step indexed by step
        (0, 1, ..., n,). Each step sub arrays is made up of the following
        variables:
    posterior : numpy array
        Array of posterior samples.
    distances : object
        Array of accepted distances.
    accepted_count : float
        Number of  posterior samples.
    trial_count : float
        Number of total samples attempted.
    epsilon : float
        Distance tolerance used.
    weights : numpy array
        Importance sampling weights. Returns an array of 1s where
        size = posterior.size when not in pmc mode.
    tau_squared : numpy array
        Gaussian kernel variances. Returns an array of 0s where
        size = posterior.size when not in pmc mode.
    eff_sample : numpy array
        Effective sample size. Returns an array of 1s where
        size = posterior.size when not in pmc mode.

    Examples
    --------
    Forth coming.

    """

    output_record = np.empty(steps, dtype=[('theta accepted', object),
                                           #('theta rejected', object),
                                           ('D accepted', object),
                                           ('n accepted', float),
                                           ('n total', float),
                                           ('epsilon', float),
                                           ('weights', object),
                                           ('tau_squared', object),
                                           ('eff sample size', object),
                                           ])

    if resume != None:
        steps = xrange(resume.size, resume.size + steps)
        output_record = stack_arrays((resume, output_record), asrecarray=True,
                                     usemask=False)
        epsilon = stats.scoreatpercentile(resume[-1]['D accepted'],
                                              per=75)
        theta = resume['theta accepted'][-1]
        weights = resume['weights'][-1]
        tau_squared = resume['tau_squared'][-1]

    else:
        steps = xrange(steps)
        epsilon = epsilon_0

    for step in steps:
        print 'Starting step {}'.format(step)
        if step == 0:
    #Fist ABC calculation
            output_record[step] = basic_abc(model, data, epsilon=epsilon,
                                            min_samples=min_samples,
                                            parallel=parallel,
                                            n_procs=n_procs, pmc_mode=False)

            theta = output_record[step]['theta accepted']
            #print theta.shape
            tau_squared = 2 * np.cov(theta)
            #print tau_squared
            weights = np.ones(theta.shape[1]) * 1.0/theta.shape[1]
            #print weights
            epsilon = stats.scoreatpercentile(output_record[step]['D accepted'],
                                              per=75)

            output_record[step]['weights'] = weights
            output_record[step]['tau_squared'] = tau_squared

            #print tau_squared
            #print weights
            #print epsilon

        else:
            #print weights, tau_squared
            #print theta
            theta_prev = theta
            weights_prev = weights
            output_record[step] = basic_abc(model, data, epsilon=epsilon,
                                            min_samples =min_samples,
                                            parallel=parallel,
                                            n_procs=n_procs, pmc_mode=True,
                                            weights=weights,
                                            theta_prev=theta_prev,
                                            tau_squared=tau_squared)

            theta = output_record[step]['theta accepted']
            epsilon = stats.scoreatpercentile(output_record[step]['D accepted'],
                                              per=75)

            #print theta


            if epsilon == 0.0:
                epsilon = 0.001

            #print theta_prev
            effective_sample = effective_sample_size(weights_prev)

            weights = calc_weights(theta_prev, theta, tau_squared, weights_prev,
                                                        prior=model.prior)

            tau_squared = 2 * weighted_covar(theta, weights)

            output_record[step]['tau_squared'] = tau_squared

            output_record[step]['eff sample size'] = effective_sample

            output_record[step]['weights'] = weights

    return output_record

#
# def calc_weights(theta_prev, theta, tau_squared, weights, prior="None"):
#
#     """
#     Calculates importance weights
#     """
#     weights_new = np.zeros_like(theta)
#
#     for i in xrange(theta.shape[0]):
#         for j in xrange(theta[i].size):
#             weights_new[i][j] = (prior[i].pdf(theta[i][j]) /
#                                 np.sum(weights[i]*stats.norm.pdf(theta[i],
#                                                                  theta_prev[i],
#                                 np.sqrt(tau_squared[0][i]))))
#
#         weights_new[i] = weights_new[i]/sum(weights_new[i])
#         #print weights_new[i]
#     return weights_new

def calc_weights(theta_prev, theta, tau_squared, weights, prior="None"):
    """
    Calculates importance weights
    """
    weights_new = np.zeros_like(weights)

    if len(theta.shape) == 1:
        norm = np.zeros_like(theta)
        for i, T in enumerate(theta):
            for j in xrange(theta_prev[0].size):
                #print T, theta_prev[0][j], tau_squared
                #print type(T), type(theta_prev), type(tau_squared)
                norm[j] = stats.norm.pdf(T, loc=theta_prev[0][j],
                                     scale=tau_squared)
            weights_new[i] = prior[0].pdf(T)/sum(weights * norm)

        return weights_new/weights_new.sum()

    else:
        norm = np.zeros(theta_prev.shape[1])
        for i in xrange(theta.shape[1]):
            prior_prob = np.zeros(theta[:, i].size)
            for j in xrange(theta[:, i].size):
                #print theta[:, i][j]
                prior_prob[j] = prior[j].pdf(theta[:, i][j])
            #assumes independent priors
            p = prior_prob.prod()

            for j in xrange(theta_prev.shape[1]):
                norm[j] = stats.multivariate_normal.pdf(theta[:, i],
                                                        mean=theta_prev[:, j],
                                                        cov=tau_squared)

            weights_new[i] = p/sum(weights * norm)

        return weights_new/weights_new.sum()

def weighted_covar(x, w):
    """

    :param x: 1 or 2 dimensional array-like, values
    :param w: 1 dimensional array-like, weights
    :return C: Weighted covariance of x or weighted variance if x is 1d
    """
    sumw = w.sum()
    assert_almost_equal(sumw, 1.0)
    if len(x.shape) == 1:
        assert x.shape[0] == w.size
    else:
        assert x.shape[1] == w.size
    sum2 = np.sum(w**2)


    if len(x.shape) == 1:
        xbar = (w*x).sum()
        var = sum(w * (x - xbar)**2)
        return var * sumw/(sumw*sumw-sum2)
    else:
        xbar = [(w*x[i]).sum() for i in xrange(x.shape[0])]
        covar = np.zeros((x.shape[0], x.shape[0]))
        for k in xrange(x.shape[0]):
            for j in xrange(x.shape[0]):
                for i in xrange(x.shape[1]):
                    covar[j,k] += (x[j,i]-xbar[j])*(x[k,i]-xbar[k]) * w[i]

        return covar * sumw/(sumw*sumw-sum2)

def effective_sample_size(w):
    '''

    :param w: array-like importance sampleing weights
    :return: float, effective sample size
    '''

    sumw = sum(w)
    sum2 = sum (w**2)
    return sumw*sumw/sum2
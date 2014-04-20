"""
Module for Approximate Bayesian Computation

Usage
-----

    ABC.py <file containing stars> <model.py> <observed>


"""
import sys
import numpy as np
import time
from scipy import stats
import pylab as plt
from simple_lib import *

#@profile



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






def main():
    start = time.time()
    N = int(sys.argv[4])

    if len(sys.argv) == 6:
        print """
        Creating and saving a realization of the forward model as {}
        """.format(sys.argv[5])

        out = file(sys.argv[5], 'w')
        N = 1


    #Load input stars and model
    model_def = __import__(sys.argv[2])

    planet_count = 0
    planet_counts, bs = [], []
    last_star = None
    for line in file(sys.argv[3]):
        line = line.split(',')
        # for k,l in enumerate(line):
        #    print k,l
        # sys.exit()
        if line[0] == 'ktc_kepler_id':
            continue
        if 1.0 >= float(line[15]) >= -1.0:
            bs.append(float(line[15]))
            planet_count += 1
            if last_star != line[0] and last_star is not None:
                planet_counts.append(planet_count)
                planet_count = 0
            last_star = line[0]
    if planet_count != 0:
        planet_counts.append(planet_count)
        #print line[15], planet_count
        #print sum(planet_counts), len(bs)

    planet_counts = np.array(planet_counts)
    bs = np.array(bs)
    #print planet_counts.sum(), bs.size
    assert planet_counts.sum() == bs.size

    stars = np.recfromcsv('stars.csv')

    ABC_model = model_def.ABC()


    epsilon = 0.01


    ac_Dc, ac_Db, ac_n, ac_sig = [], [], [], []
    re_Dc, re_Db, re_n, re_sig = [], [], [], []

    start_catalog = time.time()
    trial, accepted = 0, 0
    while accepted <= N:
        trial += 1
        #eventally this will take a vector of model parameters
        model = model_def.Model()
        star_header = model.star_parameters
        planet_header = model.planet_parameters

        #Dumb sampling for now
        binom_n = stats.randint.rvs(1, 10, 1)
        scale = stats.uniform.rvs(0, 10, 1)
        if len(sys.argv) == 6:
            binom_n = 5
            scale = 3.0  #TODO Move defaults to the model.ABC

        #Draw the number of planets per star.
        planet_numbers = model.planets_per_system(binom_n,
                                                  stars['ktc_kepler_id'].size)

        total_planets = planet_numbers.sum()

        #Initalize synthetic catalog.
        catalog = np.zeros(planet_numbers.sum(),
                           dtype={'names': star_header + planet_header,
                                  'formats': (['i8'] + ['f8'] *
                                              (len(star_header + planet_header)
                                               - 1))})

        #Draw the random model parameters.

        catalog['period'] = model.planet_period(total_planets)
        catalog['mi'] = model.mutual_inclination(scale, total_planets)
        catalog['fund_plane'] = model.fundamental_plane(total_planets)
        catalog['fund_node'] = model.fundamental_node(total_planets)
        catalog['e'] = model.eccentricity(total_planets)
        catalog['w'] = model.longitude_ascending_node(total_planets)

        for h in star_header:
            catalog[h] = np.repeat(stars[h], planet_numbers)

        assert isinstance(catalog, object)

       # print catalog.dtype.names

        #Compute derived parameters.
        catalog['a'] = semimajor_axis(catalog['period'], catalog['mass'])
        catalog['i'] = inclination(catalog['fund_plane'], catalog['mi'],
                                   catalog['fund_node'])
        catalog['b'] = impact_parameter(catalog['a'], catalog['e'],
                                        catalog['i'], catalog['w'],
                                        catalog['radius'])


        #Save a catalog
        if len(sys.argv) == 6:
            np.savetxt(out, catalog, delimiter=',',
                        fmt=['%s15']+['%10.5f']*(len(catalog.dtype.names)-1),
                        newline='\n',
                        header=','.join(x for x in star_header+planet_header),
                        comments='')
            out.close()
            sys.exit('Done making synthetic catalog.')



        #Let's do some ABC!
        transits = np.where((1.0 > catalog['b']) & (catalog['b'] > -1.0))
        pcount = np.bincount(catalog['ktc_kepler_id'][transits])

        Dc = ABC_model.distance_function(pcount[np.nonzero(pcount)],
                                         planet_counts)
        Db = ABC_model.distance_function(bs, catalog['b'][transits])

        if Dc <= epsilon and Db <= epsilon:
            accepted += 1
            ac_Dc.append(Dc)
            ac_Db.append(Db)
            ac_n.append(binom_n)
            ac_sig.append(scale)
        else:
            re_Dc.append(Dc)
            re_Db.append(Db)
            re_n.append(binom_n)
            re_sig.append(scale)

    end_time = time.time()
    print """
    {} trials of {} trials accepted
    Total: {}s  Catalog Gen + ABC: {}s
    Overhead: {}s  Time/Catalog: {}s
    """.format(N, trial, end_time - start, end_time - start_catalog,
               start_catalog - start, ((end_time - start_catalog)/float(trial)))


    f1 = plt.figure()
    plt.plot(re_n, re_sig, 'ko', alpha=.2)
    plt.plot(ac_n, ac_sig, 'o', alpha=.8, mfc='c', ms=8, mec='c')
    plt.axhline(3,ls='--')
    plt.axvline(5,ls='--')
    plt.xlabel('n', fontsize=18)
    plt.ylabel('sigma', fontsize=18)
    plt.suptitle(r'$\epsilon$ = {} n_acc = {}  N = {}'.format(epsilon,
                                                              N, trial),
                 fontsize=18)
    #plt.show()
    plt.savefig('{}.pdf'.format(sys.argv[2]))

if __name__ == "__main__":
    main()
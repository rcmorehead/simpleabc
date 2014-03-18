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
# noinspection PyNoneFunctionAssignment
def main():
    start = time.time()
    N = int(sys.argv[4])

    #Load input stars and model
    model_def = __import__(sys.argv[2])

    planet_count = 1
    planet_counts, bs = [], []
    last_star = None
    for line in file(sys.argv[3]):
        line = line.split(',')
        if line[0] == 'kepid':
            continue
        if 1.0 > float(line[7]) > -1.0:
            bs.append(float(line[7]))
            if last_star == line[0]:
                planet_count += 1
            else:
                planet_counts.append(planet_count)
                planet_count = 1
            last_star = line[0]

    planet_counts = np.array(planet_counts)
    bs = np.array(bs)

    assert planet_counts.sum() == bs.size

    stars = np.recfromcsv('stars.csv')

    model = model_def.Model()

    star_header = model.star_parameters
    planet_header = model.planet_parameters



    epsilon = 0.05


    ac_Dc, ac_Db, ac_n, ac_sig = [], [], [], []
    re_Dc, re_Db, re_n, re_sig = [], [], [], []

    start_catalog = time.time()
    for n in range(N):
        binom_n = stats.randint.rvs(1, 10, 1)
        sigma = stats.uniform.rvs(0, 10, 1)

        #Draw the number of planets per star.
        planet_numbers = model.planets_per_system(binom_n,
                                                  stars['ktc_kepler_id'].size)


        #Initalize synthetic catalog.
        catalog = np.zeros(planet_numbers.sum(),
                           dtype={'names': star_header + planet_header,
                                  'formats': (['i8'] + ['f8'] *
                                              (len(star_header + planet_header)
                                               - 1))})

        #Draw the random model parameters.
        catalog['b'] = model.planet_b(planet_numbers.sum(), sigma)
        catalog['period'] = model.planet_period(planet_numbers.sum())

        for h in star_header:
            catalog[h] = np.repeat(stars[h], planet_numbers)

        assert isinstance(catalog, object)

       # print catalog.dtype.names

        #Compute derived parameters.
        catalog['a'] = semimajor_axis(catalog['period'], catalog['mass'])

        ABC_model = model_def.ABC()

        transits = np.where((1.0 > catalog['b']) & (catalog['b'] > -1.0))
        pcount = np.bincount(catalog['ktc_kepler_id'][transits])
        Dc = ABC_model.distance_function(pcount[np.nonzero(pcount)],
                                         planet_counts)
        Db = ABC_model.distance_function(bs, catalog['b'][transits])

        if Dc <= epsilon and Db <= epsilon:
            ac_Dc.append(Dc)
            ac_Db.append(Db)
            ac_n.append(binom_n)
            ac_sig.append(sigma)
        else:
            re_Dc.append(Dc)
            re_Db.append(Db)
            re_n.append(binom_n)
            re_sig.append(sigma)

    end_time = time.time()
    print """
    Total: {}s  Catalog Gen + ABC: {}s
    Overhead: {}s  Time/Catalog: {}s
    """.format(end_time - start, end_time - start_catalog,
               start_catalog - start, ((end_time - start_catalog)/float(N)))


    f1 = plt.figure()
    plt.plot(re_n, re_sig, 'ko',alpha=.4)
    plt.plot(ac_n, ac_sig, 'bo',alpha=1)
    plt.axhline(2,ls='--')
    plt.axvline(5,ls='--')
    plt.xlabel('n', fontsize=18)
    plt.ylabel('sigma', fontsize=18)
    plt.suptitle(r'$\epsilon$ = {} n = {}'.format(epsilon,N), fontsize=18)
    #plt.show()
    plt.savefig('{}.pdf'.format(sys.argv[2]))

if __name__ == "__main__":
    main()
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

#@profile
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

    stars = np.recfromcsv('stars.csv',
                      dtype=('int','float','float','float','float'))

    model = model_def.Model()

    star_header = list(stars.dtype.names)
    planet_header = model.planet_parameters



    epsilon = 0.05


    ac_Dc, ac_Db, ac_n, ac_sig = [],[],[],[]
    re_Dc, re_Db, re_n, re_sig = [],[],[],[]

    for n in range(N):
        binom_n = stats.randint.rvs(1,10,1)
        sigma = stats.uniform.rvs(0,10,1)

        planet_numbers = model.planets_per_system(binom_n, stars['kepid'].size)
        impact_parameters = model.planet_b(planet_numbers.sum(), sigma)

        catalog = np.zeros(impact_parameters.size,
                           dtype={'names': star_header + planet_header,
                                  'formats': (['i8'] + ['f8'] * (
                                   len(star_header + planet_header) - 1))})




        catalog['b'] = impact_parameters

        for h in star_header:
            catalog[h] = np.repeat(stars[h], planet_numbers)



        assert isinstance(catalog, object)


        transits = np.where((1.0 > catalog['b']) & (catalog['b'] > -1.0))
        pcount = np.bincount(catalog['kepid'][transits])
        Dc = stats.ks_2samp(pcount[np.nonzero(pcount)], planet_counts)[0]
        Db = stats.ks_2samp(bs, catalog['b'][transits])[0]
        #print  Dc,Db





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

    print  time.time() - start
    f1 = plt.figure()
    plt.plot(re_n, re_sig, 'ko',alpha=.4)
    plt.plot(ac_n, ac_sig, 'bo',alpha=.8)
    plt.axhline(2,ls='--')
    plt.axvline(5,ls='--')
    plt.xlabel('n')
    plt.ylabel('sigma')
    plt.suptitle(r'$\epsilon$ = {} n = {}'.format(epsilon,N))
    #plt.show()
    plt.savefig('{}.eps'.format(sys.argv[2]))

if __name__ == "__main__":
    main()
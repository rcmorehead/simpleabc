"""
Module for Approximate Bayesian Computation

Usage
-----

    ABC.py <file containing stars> <model.py>


"""
import sys
import numpy as np


def main():

    #Load input stars and model
    model_def = __import__(sys.argv[2])

    f = file(sys.argv[1], 'r')
    f = f.readlines()
    star_header = f[0].strip('\n').split(',')
    stars = {head: [] for head in star_header}
    planet_header = ['rp','mp','P']
    planets = {head: [] for head in planet_header}
    model = model_def.Model()

    for i in range(1, len(f)):
        if '--' in f[i]:
            pass
        else:
            row = f[i].strip('\n').split(',')
            planet_number = model.planets_per_system()
            for n in range(planet_number):
                for h in range(len(star_header)):
                    stars[star_header[h]].append(row[h])
                for h in range(len(planet_header)):
                    planets[planet_header[h]].append(0)


    catalog = np.core.records.fromarrays([stars[x] for x in star_header] +
                                         [planets[x] for x in planet_header],
                                         names=star_header+planet_header,
                                         formats='float64, '*len(star_header
                                                                + planet_header)
                                         )

    print catalog['mass'],catalog['teff'],catalog['kepid']


if __name__ == "__main__":
    main()
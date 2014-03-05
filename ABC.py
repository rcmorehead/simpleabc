"""
Module for Approximate Bayesian Computation

Usage
-----

    ABC.py <file containing stars> <model.py> <observed>


"""
import sys
import numpy as np
import time

#@profile
def main():
    start = time.time()
    N = int(sys.argv[4])

    #Load input stars and model
    model_def = __import__(sys.argv[2])

    planet_count = 1
    planet_counts, bs = [],[]
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


    f = file(sys.argv[1], 'r')
    f = f.readlines()
    for i in range(1, len(f)):
        if '--' in f[i]:
            pass
        else:
            row = f[i].strip('\n').split(',')


    star_header = f[0].strip('\n').split(',')
    planet_header = ['rp', 'b', 'P']

    model = model_def.Model()
    planet_numbers = model.planets_per_system(binom_n)


    binom_n, sigma = 5, 2
    for i in range(N):
        stars = {head: [] for head in star_header}
        planets = {head: [] for head in planet_header}
        for i in range(1, len(f)):

            for i in range(planet_numbers[i]):
                if 1.0 > b[i] > -1.0:
                    for h in range(len(star_header)):
                        stars[star_header[h]].append(row[h])
                    for h in range(len(planet_header)):
                        if planet_header[h] == 'b':
                            planets[planet_header[h]].append(b[i])
                        else:
                            planets[planet_header[h]].append(0)
                else:
                    pass



        catalog = np.core.records.fromarrays([stars[x] for x in star_header] +
                                             [planets[x] for x in planet_header]
                                             ,names=star_header+planet_header,
                                             formats=('|S15, ' + 'float64, '
                                                        *(len(star_header
                                                        + planet_header)-1)))

        assert isinstance(catalog, object)
        print catalog['kepid'],catalog['mass'], catalog['teff'], catalog['b']

    # out = file('test_cat.dat', 'w')
    # np.savetxt(out, catalog, delimiter=',',
    #            fmt=['%s15']+['%10.5f']*(len(catalog.dtype.names)-1),
    #            newline='\n',
    #            header=','.join(x for x in star_header+planet_header),
    #            comments='')
    # out.close()

    print time.time()-start

if __name__ == "__main__":
    main()
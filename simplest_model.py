import model
from scipy import stats
import numpy as np
from simple_lib import *

class MyModel(model.Model):

        def __init__(self, stars, data):
            self.stars = stars
            self.data = data
            self.data_sum_stats = self.summary_stats(self.data)

        def draw_theta(self):
            binom_n = stats.randint.rvs(1, 10, 1)
            scale = stats.uniform.rvs(0, 10, 1)

            return binom_n, scale[0]

        def generate_data(self, theta):
            planet_numbers = (
                self.planets_per_system(theta[0],
                                        self.stars['ktc_kepler_id'].size))


            total_planets = planet_numbers.sum()

            star_header = ['ktc_kepler_id', 'teff', 'teff_err1', 'logg', 'feh',
                       'feh_err1', 'mass', 'mass_err1', 'radius', 'radius_err1',
                       'cdpp3', 'cdpp6', 'cdpp12', 'kepmag', 'days_obs']

            planet_header = ['b', 'i', 'a', 'planet_mass', 'planet_radius',
                         'period', 'mi', 'fund_plane', 'fund_node', 'e', 'w']

            #Initalize synthetic catalog.
            catalog = np.zeros(planet_numbers.sum(),
                                dtype={'names': star_header + planet_header,
                                'formats': (['i8'] + ['f8'] *
                                            (len(star_header + planet_header)
                                            - 1))})

            #Draw the random model parameters.

            catalog['period'] = self.planet_period(total_planets)
            catalog['mi'] = self.mutual_inclination(theta[1], total_planets)
            catalog['fund_plane'] = self.fundamental_plane(total_planets)
            catalog['fund_node'] = self.fundamental_node(total_planets)
            catalog['e'] = self.eccentricity(total_planets)
            catalog['w'] = self.longitude_ascending_node(total_planets)

            for h in star_header:
                catalog[h] = np.repeat(self.stars[h], planet_numbers)



            # print catalog.dtype.names

            #Compute derived parameters.
            catalog['a'] = semimajor_axis(catalog['period'], catalog['mass'])
            catalog['i'] = inclination(catalog['fund_plane'], catalog['mi'],
                                   catalog['fund_node'])
            catalog['b'] = impact_parameter(catalog['a'], catalog['e'],
                                        catalog['i'], catalog['w'],
                                        catalog['radius'])
            return catalog




        def summary_stats(self, data):
            transits = np.where((1.0 > data['b']) & (data['b'] > -1.0))
            pcount = np.bincount(data['ktc_kepler_id'][transits])
            pcount = pcount[np.where(pcount > 0)]
            bs = data['b'][transits]

            return pcount, bs

        def distance_function(self, summary_stats, summary_stats_synth):
            ksd_bi = stats.ks_2samp(summary_stats[0], summary_stats_synth[0])[0]
            ksd_sc = stats.ks_2samp(summary_stats[1], summary_stats_synth[1])[0]

            return np.sqrt(ksd_bi**2+ksd_sc**2)

        def planets_per_system(self, n, size):
            return stats.binom.rvs(n, .5, size=size)

        def planet_period(self, size):
            return 10**stats.uniform.rvs(0, 3, size=size)

        def fundamental_node(self, size):
            return stats.uniform.rvs(0, 360, size=size)

        def fundamental_plane(self, size):
            return np.degrees(np.arccos(2*stats.uniform.rvs(0, 1, size)-1))

        def mutual_inclination(self, scale, size):
            return stats.rayleigh.rvs(scale, size=size)

        def eccentricity(self, size):
            return np.zeros(size)

        def longitude_ascending_node(self, size):
            return stats.uniform.rvs(0, 360, size)



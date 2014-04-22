import model
from scipy import stats
import numpy as np
from simple_lib import *

class MyModel(model.Model):

        def draw_theta(self):
            binom_n = stats.randint.rvs(1, 10, 1)
            scale = stats.uniform.rvs(0, 10, 1)

            return binom_n, scale

        def generate_data(self, theta):
            pass

        def summary_stats(self, data):
            transits = np.where((1.0 > data['b']) & (data['b'] > -1.0))
            pcount = np.bincount(data['ktc_kepler_id'][transits])
            pcount = pcount[np.where(pcount > 0)]
            bs = data['b'][transits]

            return pcount, bs

        def distance_function(self, summary_stats, summary_stats_synth):
            return stats.ks_2samp(summary_stats, summary_stats_synth)[0]

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




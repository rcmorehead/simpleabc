import model
from scipy import stats
import numpy as np


class Model(model.PlanetsModel):

        def planets_per_system(self, n, size):
            return stats.binom.rvs(n, .5, size=size)

        def planet_b(self, n, sigma):
            return stats.norm.rvs(scale=sigma, size=n)

        def planet_period(self, size):
            return 10**stats.uniform.rvs(-2, 5, size=size)

        def fundamental_node(self, size):
            return stats.uniform.rvs(0, 360, size=size)

        def fundamental_plane(self, size):
            return np.degrees(2*stats.uniform.rvs(0, np.pi, size)-1)

        def mutual_inclination(self, scale, size):
            return stats.rayleigh.rvs(scale, size=size)


class ABC(model.ABC):

        def distance_function(self, simulation, observed):
            return stats.ks_2samp(simulation, observed)[0]

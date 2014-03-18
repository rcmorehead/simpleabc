'''
Set of base classes for generating Kepler-like synthetic catalogs for
approximate bayesian computing and various uses limited only by the user's
imagination.
'''


from scipy import stats
import numpy as np

class PlanetsModel(object):


    star_parameters = ['ktc_kepler_id', 'teff', 'teff_err1', 'logg', 'feh',
                       'feh_err1', 'mass', 'mass_err1', 'radius', 'radius_err1',
                       'cdpp3', 'cdpp6', 'cdpp12', 'kepmag', 'days_obs']
    
    planet_parameters = ['b', 'i', 'a', 'planet_mass', 'planet_radius',
                         'period', 'mi', 'fund_plane', 'fund_node']

    def __init__(self):
        pass

    def fundamental_node(self):
        pass

    def fundamental_plane(self):
        pass

    def mutual_inclination(self):
        pass

    def planet_inclination(self):
        pass

    def planet_mass(self):
        pass

    def planet_period(self):
        pass

    def planet_radius(self):
        pass

    def planets_per_system(self):
        pass


class NoiseModel(object):

    def __init__(self):
        pass


class ABC(object):

    def __init__(self):
        pass

    def abc_reduce(self, data):
        pass

    def distance_function(self, simulation, observed):
        pass



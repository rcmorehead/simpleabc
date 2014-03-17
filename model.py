'''

'''


from scipy import stats
import numpy as np

def main():
    print("model.py is not meant to be call directly")

class PlanetsModel(object):


    star_parameters = ['ktc_kepler_id', 'teff', 'teff_err1', 'logg', 'feh',
                       'feh_err1', 'mass', 'mass_err1', 'radius', 'radius_err1',
                       'cdpp3', 'cdpp6', 'cdpp12', 'kepmag', 'days_obs']
    
    planet_parameters = ['b', 'i', 'planet_mass', 'planet_radius']

    def __init__(self):
        pass

    def planets_per_system(self):
        pass

    def planet_period(self):
        pass

    def planet_mass(self):
        pass

    def planet_radius(self):
        pass

    def planet_inclination(self):
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


def main():
    pass

if __name__ == "__main__":
    main()

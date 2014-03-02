"""
Useful classes and functions for SIMPLE.
"""


class System(object):
    """
    A container class for system parameters.

    Parameters
    ----------
    #TODO
    Returns
    -------
    Nothing yet

    Methods
    -------
    System.add_planet(planet)
        Add a simple_lib.Planet object to the system

    """


    def __init__(self):
        self.planets = [] #A list may not be the best thing for this

    def add_planet(self,planet):
        self.planets.append(planet)


class Planet(object):
    """
    A container class for planet parameters.

    Probably a dead end, going to try a record array based approach

    Parameters
    ----------
    period : float
        Period in days
    mass : float
        Mass in Earth masses
    radius : float
        Radius in Earth radii
    a : float
        Semimajor axis in AU

    Attributes
    ----------
    #TODO

    Returns
    -------
    #TODO

    Methods
    -------
    planet.set_period(new_period) : float
        Set planet.period to new_period
    planet.get_period() : none
        Return planet.period
    """

    def __init__(self, period, mass, radius, a):
        self.period = period
        self.mass = mass
        self.radius = radius
        self.a = a

    def set_period(self, new_period):
        self.period = new_period

    def get_period(self):
        return self.period


def main(self):
    pass

if __name__ == "__main__":
    main()
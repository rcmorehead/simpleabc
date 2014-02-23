'''

'''

class System(object):

    def __init__(self):
        pass


class Planet(object):
    """
    A container class for planet parameters.

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

    Returns
    -------

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
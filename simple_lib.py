"""
Useful classes and functions for SIMPLE.
"""
import numpy as np


def semimajor_axis(period, mass):
    """
    Compute the semimajor axis of an object.

    This is a simple implementation of the general form Kepler's Third law.

    Parameters
    ----------
    period : int, float or numpy array
        The orbital period of the orbiting body in units of days.
    mass : int, float or array-like
        The mass of the central body (or mass sum) in units of solar mass.

    Returns
    -------
    semimajor_axis : float or numpy array
        The semimajor axis in AU.

    Examples
    --------
    >>> semimajor_axis(365.256363,1.00)
    0.999985270598628

    >>> semimajor_axis(np.linspace(1, 1000, 5),np.linspace(0.08, 4, 5))
    array([ 0.00843254,  0.7934587 ,  1.56461631,  2.33561574,  3.10657426])

    """
    return (((2.959E-4*mass)/(4*np.pi**2))*period**2.0) ** (1.0/3.0)

def main():
    pass

if __name__ == "__main__":
    main()
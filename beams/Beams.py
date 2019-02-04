import scipy as sp
import matplotlib.pyplot as plt
from scipy.special import hermite, genlaguerre
import math

class Paraxial(object):
    """Contains functions that describe different types of paraxial beams, i.e. waist, radius of curvature, etc.
    Instance variables are wavelength, waist, wavevector, and Rayleigh length (calculated from wavelength and waist)."""

    version = 1

    def __init__(self, lamb, waist):
        self.lamb = lamb
        self.k = 2 * sp.pi / lamb
        self.waist = waist
        self.RL = sp.pi * waist**2 / lamb  # Rayleigh length

    def Spot_size(self, z):
        """Returns the spot size of the beam at axial position z."""
        return self.waist * sp.sqrt(1 + (z / self.RL)**2)

    def R_curve(self, z):
        """Returns the radius of curvature of the phase front at axial position z."""
        return z + self.RL**2 / z

    def Gouy(self, z):
        """Returns the 00 mode Gouy phase at axial position z."""
        return sp.arctan(z / self.RL)

    def Trans_prof(self, x, y, z):
        """Returns the transverse profile of the beam at axial position z."""
        pass

    def ThreeD_Plotting(self, x, y, z, func):
        """Plots a 3-d surface plot of the beam's transverse profile."""
        pass

    def CM_Plotting(selfself, x, y, z, func):
        """Plots a colormap of the beam's transverse profile."""
        pass

class HG(Paraxial):

    """Creates Hermite-Gaussian beams as a TEM_n,m mode. Takes rectangular coordinates.
    This needs to be double checked for accuracy."""

    def __init__(self, lamb, waist, n, m):
        Paraxial.__init__(self, lamb, waist)
        self.n = n  # x mode
        self.m = m  # y mode

    @property
    def n(self):
        return self._n

    @n.setter
    def n(self, value):
        if isinstance(value, int)==True and value >= 0:
            self._n = value
        else:
            raise ValueError("X index must be an integer greater than zero.")
        return None

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, value):
        if isinstance(value, int)==True and value >= 0:
            self._m = value
        else:
            raise ValueError("Y index must be an integer greater than zero.")
        return None

    def Single_prof(self, nn, z, pos):
        """Returns the Hermite-Gaussian profile in a single x or y direction."""
        gouy = sp.sqrt(sp.exp(1j * (2 * nn + 1) * self.Gouy(z)) / (2**nn * math.factorial(nn) * self.Spot_size(z)))
        expon = sp.exp(-1j * self.k * z - 1j * self.k * pos**2 / (2 * self.R_curve(z)) - pos**2 / self.Spot_size(z)**2)
        return (2 / sp.pi)**(0.25) * gouy * hermite(nn)(pos) * expon

    def Trans_prof(self, x, y, z):
        print self.n
        print self._n
        return self.Single_prof(self.n, z, x) * self.Single_prof(self.m, z, y)

    def ThreeD_Plotting(self, x, y, z, func):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_surface(x, y, func(x, y, z).real, rstride=1, cstride=1, linewidth=0, cmap='ocean')
        plt.show()
        return None

    def CM_Plotting(self, x, y, z, func):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.pcolormesh(x, y, func(x, y, z).real, cmap='gray')
        plt.show()
        return None

class LG(Paraxial):

    """Create Laguerre-Gaussian beams.The radial mode is indexed by rho and the azimuthal mode is indexed by phi"""

    def __init__(self, lamb, waist, p, m):
        Paraxial.__init__(self, lamb, waist)
        self.p = p
        self.m = m

    @property
    def p(self):
        return self._p

    @p.setter
    def p(self, value):
        if isinstance(value, int) == True and value >= 0:
            self._p = value
        else:
            raise ValueError("Radial index must be an integer greater to or equal to zero")
        return None

    @property
    def m(self):
        return self._m

    @m.setter
    def m(self, value):
        if isinstance(value, int) == True and value >= 0:
            self._m = value
        else:
            raise ValueError("Azimuthal index must be an integer greater to or equal to zero")
        return None

    def Trans_prof(self, rho, phi, z):
        coef = sp.sqrt(2 * math.factorial(self._p) / ((1 + sp.kron(0, self._m)) * sp.pi * math.factorial(self._m + self._p)))
        gouy = sp.exp(1j * (2 * self._p + self._m + 1) * self.Gouy(z)) / self.Spot_size(z)
        expon = sp.exp(-1j * self.k * rho**2 / (2 * self.R_curve(z)) - rho**2 / self.Spot_size(z)**2 + 1j * self._m * phi)
        return coef * gouy * (sp.sqrt(2) * rho / self.Spot_size(z))**self._m * genlaguerre(self._p, self._m)(rho) * expon

    def ThreeD_Plotting(self, rho, theta, z, func):
        xxx, yyy = rho * sp.cos(theta), rho * sp.sin(theta)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.plot_surface(xxx, yyy, func(rho, theta, z).real, rstride=1, cstride=1, linewidth=0, cmap='ocean')
        plt.show()
        return None

    def CM_Plotting(self, rho, theta, z, func):
        xxx, yyy = rho * sp.cos(theta), rho * sp.sin(theta)
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        ax.pcolormesh(xxx, yyy, func(rho, theta, z).real, cmap='gray')
        plt.show()
        return None

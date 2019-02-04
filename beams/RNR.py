#  Reflection and Refraction
import scipy as sp

class Xwave(object):
    """Calculates reflection and transmission coefficients among other things.
    n1 (mu1) is the starting medium and n2 (mu2) is the material on which the beam is incident."""

    def __init__(self, n1, n2, mu1, mu2):
        self.n1 = n1
        self.n2 = n2
        self.mu1 = mu1
        self.mu2 = mu2
        self.crit_ang = None

    @property
    def crit_ang(self):
        return self._crit_ang

    @crit_ang.setter
    def crit_ang(self, value):
        if self.n1 <= self.n2:
            self._crit_ang = sp.pi / 2
        elif self.n1 > self.n2:
            self._crit_ang = sp.arcsin(self.n2 / self.n1)
        return None

    def Reflection(self, theta_i):
        """Calculates the reflection coefficient."""
        pass

    def Transmission(self, theta_i):
        """Calculates the transmission coefficient."""
        pass

    def Snells(self, theta_i):
        """Returns the angle of refraction given incident angel theta_i."""
        return sp.arcsin(self.n1 / self.n2 * sp.sin(theta_i))

    def Brewster(self):
        return sp.arctan(self.n2 / self.n1 * sp.sqrt((self.n1**2 - (self.mu1 / self.mu2)**2 * self.n2**2) / (self.n1**2 - self.n2**2)))

class Swave(Xwave):
    """Calculates reflection and transmission coefficients for s waves."""

    def __init__(self, n1, n2, mu1, mu2):
        Xwave.__init__(self, n1, n2, mu1, mu2)

    def Reflection(self, theta_i):
        numerator = (self.n1 * sp.cos(theta_i) - (self.mu1 / self.mu2) * sp.sqrt(self.n2**2-(self.n1 * sp.sin(theta_i))**2))
        denominator = (self.n1 * sp.cos(theta_i) + (self.mu1 / self.mu2) * sp.sqrt(self.n2**2-(self.n1 * sp.sin(theta_i))**2))
        return (numerator / denominator)**2

    def Transmission(self, theta_i):
        numerator = 2 * self.n1 * sp.cos(theta_i)
        denominator = (self.n1 * sp.cos(theta_i) + (self.mu1 / self.mu2) * sp.sqrt(self.n2**2-(self.n1 * sp.sin(theta_i))**2))
        return (numerator / denominator)**2 * sp.sqrt(1 - (self.n1 * sp.sin(theta_i) / self.n2)**2) / sp.cos(theta_i) * self.n2 / self.n1

class Pwave(Xwave):
    """Calculates reflection and transmission coefficients for p waves."""

    def __init__(self, n1, n2, mu1, mu2):
        Xwave.__init__(self, n1, n2, mu1, mu2)

    def Reflection(self, theta_i):
        numerator = (self.mu1 / self.mu2) * self.n2**2 * sp.cos(theta_i) - self.n1 * sp.sqrt(self.n2**2 - self.n1**2 * sp.sin(theta_i)**2)
        denominator = (self.mu1 / self.mu2) * self.n2**2 * sp.cos(theta_i) + self.n1 * sp.sqrt(self.n2**2 - self.n1**2 * sp.sin(theta_i)**2)
        return (numerator / denominator)**2

    def Transmission(self, theta_i):
        numerator = 2 * self.n1 * self.n2 * sp.cos(theta_i)
        denominator = (self.mu1 / self.mu2) * self.n2**2 * sp.cos(theta_i) + self.n1 * sp.sqrt(self.n2**2 - self.n1**2 * sp.sin(theta_i)**2)
        return (numerator / denominator)**2 * sp.sqrt(1 - (self.n1 * sp.sin(theta_i) / self.n2)**2) / sp.cos(theta_i) * self.n2 / self.n1

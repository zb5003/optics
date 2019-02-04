import scipy as sp
import PhysicalConstants as pc


class NormalSurface:
    """Class for finding the normal surface for uniaxial and biaxial crystals.
    For uniaxial crystals nx = ny < nz.
    For biaxial crystals nx < ny < nz."""

    def __init__(self, nx, ny, nz):
        """
        The inputted indices are sorted to comply with conventions.
        :param nx: Smallest index.
        :param ny: Second smallest index (or the same as nx if uniaxial).
        :param nz: Largest index.
        """
        self.nx = self.generate_nx([nx, ny, nz])
        self.ny = self.generate_ny([nx, ny, nz])
        self.nz = self.generate_nz([nx, ny, nz])
        self.OA = sp.arctan(self.nz / self.nx * sp.sqrt((self.ny**2 - self.nx**2) / (self.nz**2-self.ny**2)))

    def generate_nx(self, value):
        return sorted(value)[0]

    def generate_ny(self, value):
        return sorted(value)[1]

    def generate_nz(self, value):
        return sorted(value)[2]

    def n1(self, theta, phi):
        """
        Find the inner surface (slow axis) of the normal surface at an angle theta from the z-axis and phi from the x axis.
        See J. Q. Yao and T. S. Fahlen. "Calculations of optimum phase match parameters for the biaxial crystal KTiOPO4"
            Journal of Applied Physics 55, 65 (1984).
        :param theta: The polar angle, i.e. the angle between the wavevector and the z axis (nz).
        :param phi: The azimuthal angle, i.e. the angle from the x axis in the x-y plane (nx, ny).
        :return: The larger index (slow axis) of refraction for a particular propagation direction.
        """

        kx = sp.sin(theta) * sp.cos(phi)
        ky = sp.sin(theta) * sp.sin(phi)
        kz = sp.cos(theta)

        a = 1 / self.nx ** 2
        b = 1 / self.ny ** 2
        c = 1 / self.nz ** 2

        B = - kx**2 * (b + c) - ky**2 * (a + c) - kz**2 * (a + b)
        C = kx**2 * b * c + ky**2 * a * c + kz**2 * a * b

        return sp.sqrt(2) / sp.sqrt(- B + sp.sqrt(B**2 - 4 * C))

    def n2(self, theta, phi):
        """
        Find the outer surface (fast axis) of the normal surface at an angle theta from the z-axis and phi from the x axis.
        See J. Q. Yao and T. S. Fahlen. "Calculations of optimum phase match parameters for the biaxial crystal KTiOPO4"
            Journel of Applied Physics 55, 65 (1984).
        :param theta: The polar angle, i.e. the angle between the wavevector and the z axis (nz).
        :param phi: The azimuthal angle, i.e. the angle from the x axis in the x-y plane (nx, ny).
        :return: The smaller index (fast axis) of refraction for a particular propagation direction.
        """

        kx = sp.sin(theta) * sp.cos(phi)
        ky = sp.sin(theta) * sp.sin(phi)
        kz = sp.cos(theta)

        a = 1 / self.nx ** 2
        b = 1 / self.ny ** 2
        c = 1 / self.nz ** 2

        B = - kx ** 2 * (b + c) - ky ** 2 * (a + c) - kz ** 2 * (a + b)
        C = kx ** 2 * b * c + ky ** 2 * a * c + kz ** 2 * a * b

        return sp.sqrt(2) / sp.sqrt(- B - sp.sqrt(B ** 2 - 4 * C))

    def eigen_polarization(self, theta, phi, n):
        """
        Find the eigenpolarization for a wave propagating in a particular direction in an anisotropic medium.
        :param theta: The polar angle, i.e. the angle between the wavevector and the z axis (nz).
        :param phi: The azimuthal angle, i.e. the angle from the x axis in the x-y plane (nx, ny).
        :param n: The index of refraction for the polarization and direction of propagation.
        :return: An array containing the unit vector pointing in the polarization direction.
        """
        x = sp.sin(theta) * sp.cos(phi) / (n**2 - self.nx**2)
        y = sp.sin(theta) * sp.sin(phi) / (n**2 - self.ny**2)
        z = sp.cos(theta) / (n**2 - self.ny**2)

        return sp.asarray([x, y, z])


class PhaseMatch:

    def __init__(self, nfund_x, nfund_y, nfund_z, nshg_x, nshg_y, nshg_z):
        """
        Parent class for finding the SHG phase matching values of phi and theta.
        J. Q. Yao and T. S. Fahlen. "Calculations of optimum phase match parameters for the biaxial crystal KTiOPO4"
            Journal of Applied Physics 55, 65 (1984)
        M. V. Hoben. "Phase-Matched Second-Harmonic Generation in Biaxial Crystals." Journal of Applied Physics 38, 4365
            (1967).
        :param nfund_x: nx for the fundamental frequency.
        :param nfund_y: ny for the fundamental frequency.
        :param nfund_z: nz for the fundamental frequency.
        :param nshg_x: nx for the second harmonic.
        :param nshg_y: ny for the second harmonic.
        :param nshg_z: nz for the second harmonic.
        """
        self.fund = NormalSurface(nfund_x, nfund_y, nfund_z)
        self.shg = NormalSurface(nshg_x, nshg_y, nshg_z)
        self.error = 0.0000001
        self.theta_limit = sp.pi / 2**34
        self.nphi = 2000

    def phase_match_curve(self):
        """
        Return an array of phi's and theta's such that the phase matching condition is met.
        :return:
        """
        pass

    def from_equator(self, phi):
        """Algorithm for finding the phase matching condition in a particular slice of phi starting at theta = pi / 2."""
        pass

    def from_pole(self, phi):
        """Algorithm for finding the phase matching condition in a particular slice of phi starting at theta = pi / 2."""
        pass


class TypeI(PhaseMatch):

    def __init__(self, nfund_x, nfund_y, nfund_z, nshg_x, nshg_y, nshg_z):
        """
        Find the phase matching parameters phi and theta for Type I phase matching, nshg_1 = nfund_2.
        :param nfund_x: nx for the fundamental frequency.
        :param nfund_y: ny for the fundamental frequency.
        :param nfund_z: nz for the fundamental frequency.
        :param nshg_x: nx for the second harmonic.
        :param nshg_y: ny for the second harmonic.
        :param nshg_z: nz for the second harmonic.
        """
        PhaseMatch.__init__(self, nfund_x, nfund_y, nfund_z, nshg_x, nshg_y, nshg_z)

    def phase_match_curve(self):
        """
        Finds the polar angles that satisfy the type I phase matching conditions for azimuthal angles 0 to pi/2.  For a
                given azimuthal angle, there might be up to two phase matching polar angles.
        :return: 3xN array containing the phase matching angles. The first row are the azimuthal angle, the second row
                are the larger polar angles (if they exist), and the third row are the smaller polar angles (if they
                exist).
        """
        phi = sp.linspace(0, sp.pi / 2, self.nphi, endpoint=True)
        parameters = sp.zeros((3, self.nphi))

        for i in range(self.nphi):
            parameters[0, i] = phi[i]
            parameters[1, i] = self.from_equator(phi[i])
            parameters[2, i] = self.from_pole(phi[i])

        return parameters

    def from_equator(self, phi):
        """
        A search algorithm to find the polar angle that produces type I phase matching at a given azimuthal angle.
        The algorithm starts at theta=pi/2 (the equator).
        :param phi: Azimuthal angle at which the phase matching polar angle is desired.
        :return: The polar angle that, together with phi, will satisfy the phase matching conditions. If none exists,
                 then nan is returned.
        """
        theta1 = sp.pi / 2
        nshg1 = self.shg.n1(theta1, phi)
        nfund1 = self.fund.n2(theta1, phi)
        delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))

        nshg2 = self.shg.n1((self.shg.OA + self.fund.OA) / 2, phi).real
        nfund2 = self.fund.n2((self.shg.OA + self.fund.OA) / 2, phi)
        delta_2 = int((nshg2 - nfund2) / abs(nshg2 - nfund2))
        if (delta_1 == delta_2):
            theta1 = sp.nan

        else:
            direction = -1
            error = abs(nshg1 - nfund1)
            theta_step = sp.pi / 2 - self.shg.OA
            while (error > self.error):
                delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                theta_temp = theta1
                theta1 = theta1 + direction * theta_step / 2
                theta_step = theta_step / 2
                nshg1 = self.shg.n1(theta1, phi)
                nfund1 = self.fund.n2(theta1, phi)
                delta_2 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                if (delta_1 != delta_2):
                    direction = int((theta_temp - theta1) / abs(theta_temp - theta1))
                else:
                    pass
                error = abs(nshg1 - nfund1)

        return theta1

    def from_pole(self, phi):
        """
        A search algorithm to find the polar angle that produces type I phase matching at a given azimuthal angle.
        The algorithm starts at theta=0 (the pole).
        :param phi: Azimuthal angle at which the phase matching polar angle is desired.
        :return: The polar angle that, together with phi, will satisfy the phase matching conditions. If none exists,
                 then nan is returned.
        """
        theta1 = 0
        nshg1 = self.shg.n1(theta1, phi)
        nfund1 = self.fund.n2(theta1, phi)
        delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))

        nshg2 =self.shg.n1((self.shg.OA + self.fund.OA) / 2, phi).real
        nfund2 = self.fund.n2((self.shg.OA + self.fund.OA) / 2, phi)
        delta_2 = int((nshg2 - nfund2) / abs(nshg2 - nfund2))
        if (delta_1 == delta_2):
            theta1 = sp.nan

        else:
            direction = +1
            error = abs(nshg1 - nfund1)
            theta_step = self.shg.OA + 0.1
            while (error > self.error):
                delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                theta_temp = theta1
                theta1 = theta1 + direction * theta_step / 2
                theta_step = theta_step / 2
                nshg1 = self.shg.n1(theta1, phi)
                nfund1 = self.fund.n2(theta1, phi)
                delta_2 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                if (delta_1 != delta_2):
                    direction = int((theta_temp - theta1) / abs(theta_temp - theta1))
                else:
                    pass
                error = abs(nshg1 - nfund1)

        return theta1


class TypeII(PhaseMatch):

    def __init__(self, nfund_x, nfund_y, nfund_z, nshg_x, nshg_y, nshg_z):
        """
        Find the phase matching parameters phi and theta for Type II phase matching, nshg_1 = 1 / 2 (nfund_2 + nfund_1).
        :param nfund_x: nx for the fundamental frequency.
        :param nfund_y: ny for the fundamental frequency.
        :param nfund_z: nz for the fundamental frequency.
        :param nshg_x: nx for the second harmonic.
        :param nshg_y: ny for the second harmonic.
        :param nshg_z: nz for the second harmonic.
        """
        Phase_Match.__init__(self, nfund_x, nfund_y, nfund_z, nshg_x, nshg_y, nshg_z)

    def phase_match_curve(self):
        """
        Finds the polar angles that satisfy the type II phase matching conditions for azimuthal angles 0 to pi/2.  For a given
                  azimuthal angle, there might be up to two phase matching polar angles.
        :return: 3xN array containing the phase matching angles. The first row are the azimuthal angle, the second row
                are the larger polar angles (if they exist), and the third row are the smaller polar angles (if they
                exist).
        """
        phi = sp.linspace(0, sp.pi / 2, self.nphi, endpoint=True)
        parameters = sp.zeros((3, self.nphi))

        for i in range(self.nphi):
            parameters[0, i] = phi[i]
            parameters[1, i] = self.from_equator(phi[i])
            parameters[2, i] = self.from_pole(phi[i])

        return parameters

    def from_equator(self, phi):
        """
        A search algorithm to find the polar angle that produces type II phase matching at a given azimuthal angle.
        The algorithm starts at theta=pi/2 (the equator).
        :param phi: Azimuthal angle at which the phase matching polar angle is desired.
        :return: The polar angle that, together with phi, will satisfy the phase matching conditions. If none exists,
                 then nan is returned.
        """
        theta1 = sp.pi / 2
        nshg1 = self.shg.n1(theta1, phi)
        nfund1 = 0.5 * (self.fund.n2(theta1, phi) + self.fund.n1(theta1, phi))
        delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))

        nshg2 = self.shg.n1((self.shg.OA + self.fund.OA) / 2, phi).real
        nfund2 = 0.5 * (self.fund.n1((self.shg.OA + self.fund.OA) / 2, phi) + self.fund.n2((self.shg.OA + self.fund.OA) / 2, phi))
        delta_2 = int((nshg2 - nfund2) / abs(nshg2 - nfund2))
        if (delta_1 == delta_2):
            theta1 = sp.nan

        else:
            direction = -1
            error = abs(nshg1 - nfund1)
            theta_step = sp.pi / 2 - self.shg.OA
            while (error > self.error):
                delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                theta_temp = theta1
                theta1 = theta1 + direction * theta_step / 2
                theta_step = theta_step / 2
                nshg1 = self.shg.n1(theta1, phi)
                nfund1 = 0.5 * (self.fund.n2(theta1, phi) + self.fund.n1(theta1, phi))
                delta_2 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                if (delta_1 != delta_2):
                    direction = int((theta_temp - theta1) / abs(theta_temp - theta1))
                else:
                    pass
                error = abs(nshg1 - nfund1)

        return theta1

    def from_pole(self, phi):
        """
        A search algorithm to find the polar angle that produces type II phase matching at a given azimuthal angle.
        The algorithm starts at theta=0 (the pole).
        :param phi: Azimuthal angle at which the phase matching polar angle is desired.
        :return: The polar angle that, together with phi, will satisfy the phase matching conditions. If none exists,
                 then nan is returned.
        """
        theta1 = 0
        nshg1 = self.shg.n1(theta1, phi)
        nfund1 = 0.5 * (self.fund.n1(theta1, phi) + self.fund.n2(theta1, phi))
        delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))

        nshg2 = self.shg.n1((self.shg.OA + self.fund.OA) / 2, phi).real
        nfund2 = 0.5 * (self.fund.n1((self.shg.OA + self.fund.OA) / 2, phi) + self.fund.n2((self.shg.OA + self.fund.OA) / 2, phi))
        delta_2 = int((nshg2 - nfund2) / abs(nshg2 - nfund2))
        if (delta_1 == delta_2):
            theta1 = sp.nan

        else:
            direction = +1
            error = abs(nshg1 - nfund1)
            theta_step = self.shg.OA + 0.1
            while (error > self.error):
                delta_1 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                theta_temp = theta1
                theta1 = theta1 + direction * theta_step / 2
                theta_step = theta_step / 2
                nshg1 = self.shg.n1(theta1, phi)
                nfund1 = 0.5 * (self.fund.n1(theta1, phi) + self.fund.n2(theta1, phi))
                delta_2 = int((nshg1 - nfund1) / abs(nshg1 - nfund1))
                if (delta_1 != delta_2):
                    direction = int((theta_temp - theta1) / abs(theta_temp - theta1))
                else:
                    pass
                error = abs(nshg1 - nfund1)

        return theta1






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


class SFG:
    """To add: 1.) some general form of shgpmI_theta/shg_pmII_theta (i.e. if multivalued or for sfg).
               2.) method to calculate deff.
               3.) note that lambd must be in microns"""

    def sellmeier(self):
        """Generic form of Sellmeier equation for the material."""
        raise NotImplementedError

    def sellmeier_x(self, lambd):
        """
        Calculate nx based on te Sellmeier equation.
        :param lambd: Wavelength in microns.
        :return: Index
        """
        raise NotImplementedError

    def sellmeier_y(self, lambd):
        """
        Calculate ny based on te Sellmeier equation.
        :param lambd: Wavelength in microns.
        :return: Index
        """
        raise NotImplementedError

    def sellmeier_z(self, lambd):
        """
        Calculate nx based on te Sellmeier equation.
        :param lambd: Wavelength in microns.
        :return: Index
        """
        raise NotImplementedError

    def inner_normal(self, lambd, theta, phi):
        """
        Return the index for a material for the fast axis (inner surface of normal surface) for a given propagation direction.
        This method uses the class Normal_Surface.
        :param lambd: Wavelength in microns.
        :param theta: Polar angle (from the z axis).
        :param phi: Azimuthal angle (in the x-y plane).
        :return: The index of refraction given by the inner surface (slow axis) of the normal surface for a given
                 propagation direction.
        """
        nx = self.sellmeier_x(lambd)
        ny = self.sellmeier_y(lambd)
        nz = self.sellmeier_z(lambd)

        norm1 = NormalSurface(nx, ny, nz)

        return norm1.n1(theta, phi)

    def outer_normal(self, lambd, theta, phi):
        """
        Return the index for a material for the slow axis (outer surface of normal surface) for a given propagation direction.
        This method uses the class Normal_Surface.
        :param lambd: Wavelength in microns.
        :param theta: Polar angle (from the z axis).
        :param phi: Azimuthal angle (in the x-y plane).
        :return: The index of refraction given by the outer surface (fast axis) of the normal surface for a given
                 propagation direction.
        """
        nx = self.sellmeier_x(lambd)
        ny = self.sellmeier_y(lambd)
        nz = self.sellmeier_z(lambd)

        norm1 = NormalSurface(nx, ny, nz)

        return norm1.n2(theta, phi)

    def gaussian_pulse(self, center, width, amplitude, lambd):
        """
        The amplitude of a frequency component c/lambd in a gaussian pulse of width width amplitude amplitude and
        centered at center.
        :param center: Central wavelength.
        :param width: Pulse width (in wavelength).
        :param amplitude: Pulse amplitude.
        :param lambd: Wavelength at which the amplitude is sought.
        :return: Amplitude of the component with wavelength lambd.
        """
        return amplitude * sp.exp(-(lambd - center)**2 / (2 * width**2)) / sp.sqrt(2 * width**2 * sp.pi)

    def delta_k_II(self, lambd1, lambd3, theta, phi):
        """
        Calculate the phase mismatch, k1 + k2 - k3 for lambd1 < lambd3 and for type II phase matching and a given
        propagation direction.
        :param lambd1: One of the fundamental wavelengths.
        :param lambd3: The target sum frequency.
        :param theta: Polar angle.
        :param phi: Azimuthal angle.
        :return: The absolute value of the phase mismatch.
        """
        lambd2 = lambd1 * lambd3 / (lambd1 - lambd3)  # From w3 = w1 + w2, w=2pi*c/lambda
        n1 = self.inner_normal(lambd1, theta, phi)
        n2 = self.outer_normal(lambd2, theta, phi)
        n3 = self.inner_normal(lambd3, theta, phi)

        return abs(2 * sp.pi * ((n3 - n2) / lambd3 + (n2 - n1) / lambd1))

    def sum_frequenc_II(self, I1, I2, lambd1, lambd3, L, theta, phi):
        """
        Calculate the intensity of the summed frequency pulse with fundamental pulse intensities I1 and I2 for Type II
        phase matching.
        Equation 2.2.19 on page 78 of Nonlinear Optics by Robert Boyd, Third Edition.
        :param I1: Intensity of one of the fundamental beams.
        :param I2: Intensity of the other fundamental beam.
        :param lambd1: One of the fundamental wavelengths.
        :param lambd3: The target sum frequency.
        :param L: The length of the crystal.
        :param theta: Polar angle.
        :param phi: Azimuthal angle.
        :return: Intensity of the type II summed frequency beam.
        """
        omega3 = pc.c / lambd3
        lambd2 = lambd1 * lambd3 / (lambd1 - lambd3)  # From w3 = w1 + w2, w=2pi*c/lambda

        deff = 7.36e-12

        n1 = self.outer_normal(lambd1, theta, phi)
        n2 = self.inner_normal(lambd2, theta, phi)
        n3 = self.inner_normal(lambd3, theta, phi)

        del_k = self.delta_k_II(lambd1, lambd3, theta, phi)

        return 8 * deff ** 2 * omega3 ** 2 * I1 * I2 / (n1 * n2 * n3 * pc.epsilon0 * pc.c ** 2) * L ** 2 * sp.sinc(del_k * L / 2) ** 2

    def delta_k_I(self, lambd1, lambd3, theta, phi):
        """
        Calculate the phase mismatch, k1 + k2 - k3 for lambd1 < lambd3 and for type I phase matching and a given
        propagation direction.
        :param lambd1: One of the fundamental wavelengths.
        :param lambd3: The target sum frequency.
        :param theta: Polar angle.
        :param phi: Azimuthal angle.
        :return: The absolute value of the phase mismatch.
        """
        lambd2 = lambd1 * lambd3 / (lambd1 - lambd3)  # From w3 = w1 + w2, w=2pi*c/lambda
        n1 = self.outer_normal(lambd1, theta, phi)
        n2 = self.outer_normal(lambd2, theta, phi)
        n3 = self.inner_normal(lambd3, theta, phi)

        return abs(2 * sp.pi * ((n3 - n2) / lambd3 + (n2 - n1) / lambd1))

    def sum_frequenc_I(self, I1, I2, lambd1, lambd3, L, theta, phi):
        """
        Calculate the intensity of the summed frequency pulse with fundamental pulse intensities I1 and I2 for Type I
        phase matching.
        Equation 2.2.19 on page 78 of Nonlinear Optics by Robert Boyd, Third Edition.
        :param I1: Intensity of one of the fundamental beams.
        :param I2: Intensity of the other fundamental beam.
        :param lambd1: One of the fundamental wavelengths.
        :param lambd3: The target sum frequency.
        :param L: The length of the crystal.
        :param theta: Polar angle.
        :param phi: Azimuthal angle.
        :return: Intensity of the type I summed frequency beam.
        """
        omega3 = pc.c / lambd3
        lambd2 = lambd1 * lambd3 / (lambd1 - lambd3)  # From w3 = w1 + w2, w=2pi*c/lambda

        deff = 7.36e-12  # this is really the type II deff

        n1 = self.outer_normal(lambd1, theta, phi)
        n2 = self.outer_normal(lambd2, theta, phi)
        n3 = self.inner_normal(lambd3, theta, phi)

        del_k = self.delta_k_I(lambd1, lambd3, theta, phi)

        return 8 * deff**2 * omega3**2 * I1 * I2 / (n1 * n2 * n3 * pc.epsilon0 * pc.c**2) * L**2 * sp.sinc(del_k * L / 2)**2

    def shg_perfect_type_I(self, I1, lambd1, L, theta, phi):
        """
        Calculate the amplitude of the sh with perfect type I phase matching.  Effectively, there is only one incident
        beam since both fundamental beams have the same intensity and wavelength.
        Equation 2.2.19 on page 78 of Nonlinear Optics by Robert Boyd, Third Edition.
        :param I1: Intensity of the fundamental.
        :param lambd1: The wavelength of the fundamental.
        :param L: Length of the crystal.
        :param theta: Polar angle.
        :param phi: Azimuthal angle.
        :return: Intensity of the type I second harmonic.
        """

        deff = 7.36e-12

        n1 = self.outer_normal(lambd1, theta, phi)
        n2 = self.outer_normal(lambd1, theta, phi)
        n3 = self.inner_normal(2 * lambd1, theta, phi)

        omega3 = pc.c * n1 / lambd1 + pc.c * n2 / lambd1

        del_k = self.delta_k_I(lambd1, lambd1 / 2, theta, phi)

        return 8 * deff ** 2 * omega3 ** 2 * I1**2 / (n1 * n2 * n3 * pc.epsilon0 * pc.c**2) * L**2 * sp.sinc(del_k * L / 2)**2



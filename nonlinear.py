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
class KTP(SFG):
    """Class containing the material properties of Potassium Titanyl Phosphate (KTP)."""

    coef = {'Ax': 3.0129, 'Bx': 0.03807, 'Cx': 0.04283, 'Dx': 0.01664,
            'Ay': 3.0333, 'By': 0.04106, 'Cy': 0.04946, 'Dy': 0.01659,
            'Az': 3.3209, 'Bz': 0.05305, 'Cz': 0.05960, 'Dz': 0.01763}

    def __init__(self):
        SFG.__init__(self)

    def sellmeier(self, lambd, A, B, C, D):
        # lambd = lambd * 10**6
        return sp.sqrt(A + B / (lambd ** 2 - C) - D * lambd ** 2)

    def sellmeier_x(self, lambd):
        return self.sellmeier(lambd, KTP.coef['Ax'], KTP.coef['Bx'], KTP.coef['Cx'], KTP.coef['Dx'])

    def sellmeier_y(self, lambd):
        return self.sellmeier(lambd, KTP.coef['Ay'], KTP.coef['By'], KTP.coef['Cy'], KTP.coef['Dy'])

    def sellmeier_z(self, lambd):
        return self.sellmeier(lambd, KTP.coef['Az'], KTP.coef['Bz'], KTP.coef['Cz'], KTP.coef['Dz'])

    def shg_pmI_theta(self, phi, lambd_fund, lambd_sh):
        """Returns the polar angle required for type I phase matching for a given azimuthal angle."""
        nx_fund = self.sellmeier_x(lambd_fund)
        ny_fund = self.sellmeier_y(lambd_fund)
        nz_fund = self.sellmeier_z(lambd_fund)

        nx_sh = self.sellmeier_x(lambd_sh)
        ny_sh = self.sellmeier_y(lambd_sh)
        nz_sh = self.sellmeier_z(lambd_sh)

        typeI = TypeI(nx_fund, ny_fund, nz_fund, nx_sh, ny_sh, nz_sh)

        return typeI.from_equator(phi)

    def shg_pmII_theta(self, phi, lambd_fund, lambd_sh):
        """Returns the polar angle required for type II phase matching for a given azimuthal angle."""
        nx_fund = self.sellmeier_x(lambd_fund)
        ny_fund = self.sellmeier_y(lambd_fund)
        nz_fund = self.sellmeier_z(lambd_fund)

        nx_sh = self.sellmeier_x(lambd_sh)
        ny_sh = self.sellmeier_y(lambd_sh)
        nz_sh = self.sellmeier_z(lambd_sh)

        typeII = TypeII(nx_fund, ny_fund, nz_fund, nx_sh, ny_sh, nz_sh)

        return typeII.from_equator(phi)


class YSO(SFG):
    """Class containing the material properties of Potassium Titanyl Phosphate (KTP).  The Sellmeier coefficients found
    in YSO.coef are valid in for the wavelength range0.4358-0.6438 microns and can be found in:

    R. Beach, M. D. Shinn, L. Davis, R. W. Solarz, and W. F. Krupke. "Optical Absorption and Stimulated Emission of
        Neodymium in Yttrium Orthosilicate IEEE J. Quantum Electron 26, 1405 (1990)."""

    coef = {'Ax': 3.0895, 'Bx': 0.0334, 'Cx': 0.0043, 'Dx': 0.0199,
            'Ay': 3.1173, 'By': 0.0283, 'Cy': -0.0133, 'Dy': 0.00,
            'Az': 3.1871, 'Bz': 0.03022, 'Cz': -0.0138, 'Dz': 0.00}

    def __init__(self):
        SFG.__init__(self)

    def sellmeier(self, lambd, A, B, C, D):
        return sp.sqrt(A + B / (lambd ** 2 + C) + D * lambd ** 2)

    def sellmeier_x(self, lambd):
        return self.sellmeier(lambd, YSO.coef['Ax'], YSO.coef['Bx'], YSO.coef['Cx'], YSO.coef['Dx'])

    def sellmeier_y(self, lambd):
        return self.sellmeier(lambd, YSO.coef['Ay'], YSO.coef['By'], YSO.coef['Cy'], YSO.coef['Dy'])

    def sellmeier_z(self, lambd):
        return self.sellmeier(lambd, YSO.coef['Az'], YSO.coef['Bz'], YSO.coef['Cz'], YSO.coef['Dz'])

    def shg_pmI_theta(self, phi, lambd_fund, lambd_sh):
        """Returns the polar angle required for type I phase matching for a given azimuthal angle."""
        nx_fund = self.sellmeier_x(lambd_fund)
        ny_fund = self.sellmeier_y(lambd_fund)
        nz_fund = self.sellmeier_z(lambd_fund)

        nx_sh = self.sellmeier_x(lambd_sh)
        ny_sh = self.sellmeier_y(lambd_sh)
        nz_sh = self.sellmeier_z(lambd_sh)

        typeI = TypeI(nx_fund, ny_fund, nz_fund, nx_sh, ny_sh, nz_sh)

        return typeI.from_equator(phi)

    def shg_pmII_theta(self, phi, lambd_fund, lambd_sh):
        """Returns the polar angle required for type II phase matching for a given azimuthal angle."""
        nx_fund = self.sellmeier_x(lambd_fund)
        ny_fund = self.sellmeier_y(lambd_fund)
        nz_fund = self.sellmeier_z(lambd_fund)

        nx_sh = self.sellmeier_x(lambd_sh)
        ny_sh = self.sellmeier_y(lambd_sh)
        nz_sh = self.sellmeier_z(lambd_sh)

        typeII = TypeII(nx_fund, ny_fund, nz_fund, nx_sh, ny_sh, nz_sh)

        return typeII.from_equator(phi)


class Sellmeier:
    """
    Use a Sellmeier equation of the form

            n**2 = A + B / (lambda**2 + C) + D * lambda**2

    to calculate the various indicies of refraction for a material.
    Lambda should be in microns.
    """

    def __init__(self, coef):
        """
        Initialize the an instance.
        :param coef: Dictionary containing the Sellmeier coefficients for each axis.
        """
        self.coef = coef

    def nx(self, lambd):
        """
        Calculate the index of refraction for the x-axis of the material for a given wavelength (in microns)
        :param lambd: The wavelength in microns.
        :return: The index of refraction along the x-axis of the material.
        """
        return sp.sqrt(
            self.coef['Ax'] + self.coef['Bx'] / (lambd ** 2 + self.coef['Cx']) + self.coef['Dx'] * lambd ** 2)

    def ny(self, lambd):
        """
        Calculate the index of refraction for the y-axis of the material for a given wavelength (in microns)
        :param lambd: The wavelength in microns.
        :return: The index of refraction along the y-axis of the material.
        """
        return sp.sqrt(
            self.coef['Ay'] + self.coef['By'] / (lambd ** 2 + self.coef['Cy']) + self.coef['Dy'] * lambd ** 2)

    def nz(self, lambd):
        """
        Calculate the index of refraction for the z-axis of the material for a given wavelength (in microns)
        :param lambd: The wavelength in microns.
        :return: The index of refraction along the z-axis of the material.
        """
        return sp.sqrt(
            self.coef['Az'] + self.coef['Bz'] / (lambd ** 2 + self.coef['Cz']) + self.coef['Dz'] * lambd ** 2)

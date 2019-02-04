#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   Program for calculating and plotting the coherence length for SFG adding to a particular wavelength when the
# incident angle corresponds to the phase matching conditions for SHG.  The point of this program was to analyze the
# effects of bandwidth when doubling an ultrashort pulse.  This program is only applicable to type II phase matching.
# The program proceeds as such:
#
# 1.) Definition of physical parameters.
# 2.) Initialize an instance of KTP (or whatever crystal you have).
# 3.) Find the principle indices of refracion for the center frequency of the fundamental pulse and the second harmonic.
# 4.) Find the polar angle appropriate for type II phase matching for a specified azimuthal angle.
# 5.) Find the fast index (lower index) for the target second harmonic for propagation in the (second harmonic) phase
#     matching direction.
# 6.) Find the indices of refraction for the two fundamental wavelengths and calculate the coherence length.
# 7.) Plot.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from Birefringence import *
import matplotlib.pyplot as plt
import scipy as sp

# 1.)___________________________________________________________________________________________________________________
points = 200  # number of points to plot
lambd3 = 0.532  # target wavelength
lambd1_center = 1.064  # center of the ultrashort pulse
lambd1 = sp.linspace(1.040, 1.100, points)  # wavelength components in the ultrashort pulse
lambd2 = lambd1 * lambd3 / (lambd1 - lambd3)  # co-components of the ultrashort pulse that, when added to lambd1, will
                                              # give lambd3
L = sp.zeros(points)  # array to hold coherence lengths

# 2.)___________________________________________________________________________________________________________________
get_nxyz = KTP()

# 3.)___________________________________________________________________________________________________________________
n3_x = get_nxyz.sellmeier_x(lambd3)
n3_y = get_nxyz.sellmeier_y(lambd3)
n3_z = get_nxyz.sellmeier_z(lambd3)
n1_center_x = get_nxyz.sellmeier_x(lambd1_center)
n1_center_y = get_nxyz.sellmeier_y(lambd1_center)
n1_center_z = get_nxyz.sellmeier_z(lambd1_center)

# 4.)___________________________________________________________________________________________________________________
KTPII = TypeII(n1_center_x, n1_center_y, n1_center_z, n3_x, n3_y, n3_z)
phi = 30 * sp.pi / 180
theta = KTPII.from_equator(phi)

# 5.)___________________________________________________________________________________________________________________
n3_normal = Normal_Surface(n3_x, n3_y, n3_z)
n3 = n3_normal.n1(theta, phi)  # The higher frequency is always the fast polarization

# 6.)___________________________________________________________________________________________________________________
for i in range(points):
    n1_x = get_nxyz.sellmeier_x(lambd1[i])
    n1_y = get_nxyz.sellmeier_y(lambd1[i])
    n1_z = get_nxyz.sellmeier_z(lambd1[i])
    n2_x = get_nxyz.sellmeier_x(lambd2[i])
    n2_y = get_nxyz.sellmeier_y(lambd2[i])
    n2_z = get_nxyz.sellmeier_z(lambd2[i])

    n1_normal = Normal_Surface(n1_x, n1_y, n1_z)
    n1 = n1_normal.n1(theta, phi)

    n2_normal = Normal_Surface(n2_x, n2_y, n2_z)
    n2 = n2_normal.n2(theta, phi)

    L[i] = lambd1[i] * lambd3 / (2 * ((n3 - n2) * lambd1[i] + (n2 - n1) * lambd3))

    del n1_normal
    del n2_normal

# 7.)___________________________________________________________________________________________________________________
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(lambd1, abs(L) / 1000)
ax1.axhline(4, color='red', label="Crystal Length")
ax1.axvline(1.064, color='green', label=r"Center $\lambda_1$ (1.064 $\mu m$)")
ax1.set_xlim([1.05, 1.08])
ax1.set_ylim([0, 10])
ax1.set_title("Coherence Length for SFG Within Linewidth")
ax1.set_xlabel(r"$\lambda_1$ ($\mu m$)")
ax1.set_ylabel(r"Coherence Length ($mm$)")
ax1.legend(bbox_to_anchor=(1, 1), loc=2)
# plt.savefig("Coherence Length as the Fundamental Frequency is Detuned for phi=30.pdf", bbox_inches="tight")
# plt.close()
plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   This program finds the SHG phase matching angles for type I and II phase matching for the first octant (0 < theta <
# pi/2 and 0 < phi < pi/2 where theta is the polar angle and phi is the azimuthal angle).  It then plots these angles.
# It proceeds as follows:
#
# 1.) Define physical parameters for the problem.
# 2.) Initialize an instance of the YSO class (or whatever crystal class).
# 3.) Find the principle indices of refraction for the fundamental and second harmonic.
# 4.) Initialize instances of classes TypeI and TypeII (for calculating the phase matching angles).
# 5.) Find the phase matching angles in the first octant.
# 6.) Plot.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from Birefringence import *
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np

# 1.)___________________________________________________________________________________________________________________
lambd1 = 1.055  # fundamental wavelength in microns
lambd3 = 0.527  # second harmonic wavelength in microns
angle = sp.linspace(0, sp.pi / 2, 500)

# 2.)___________________________________________________________________________________________________________________
n = YSO()

# 3.)___________________________________________________________________________________________________________________
nx_fund = n.sellmeier_x(lambd1)
ny_fund = n.sellmeier_y(lambd1)
nz_fund = n.sellmeier_z(lambd1)
nx_harm = n.sellmeier_x(lambd3)
ny_harm = n.sellmeier_y(lambd3)
nz_harm = n.sellmeier_z(lambd3)

# print(np.round(n.sellmeier_x(0.527), 3), n.sellmeier_y(0.527), n.sellmeier_z(0.527))
# print(n.sellmeier_x(1.055), n.sellmeier_y(1.055), n.sellmeier_z(1.055))
# print(n.sellmeier_x(0.6328), n.sellmeier_y(0.6328), n.sellmeier_z(0.6328))
# print(n.shg_pmI_theta(sp.pi/2, 1.064, 0.532)*180 / sp.pi, n.shg_pmII_theta(sp.pi/2, 1.064, 0.532)*180 / sp.pi)

# 4.)___________________________________________________________________________________________________________________
KTPI = TypeI(nx_fund, ny_fund, nz_fund, nx_harm, ny_harm, nz_harm)
KTPII = TypeII(nx_fund, ny_fund, nz_fund, nx_harm, ny_harm, nz_harm)

# 5.)___________________________________________________________________________________________________________________
matchedI = KTPI.phase_match_curve() * 180 / 3.14159
matchedII = KTPII.phase_match_curve() * 180 / 3.14159

# 6.)___________________________________________________________________________________________________________________
fig = plt.figure(figsize=(15, 11))
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(matchedI[0, :], matchedI[1, :], label="Type I phase matching", color="blue")
ax1.plot(matchedI[0, :], matchedI[2, :])
ax1 = fig.add_subplot(1, 1, 1)
ax1.plot(matchedII[0, :], matchedII[1, :], label="Type II phase matching", color="green")
ax1.plot(matchedII[0, :], matchedII[2, :])
ax1.text(95.5, 88, r"527 nm: n$_x$="+str(np.round(n.sellmeier_x(lambd3), 4))+"; n$_y$="
         +str(np.round(n.sellmeier_y(lambd3), 4))+"; n$_z$="+str(np.round(n.sellmeier_z(lambd3), 4))+"\n \n"
         r"1055 nm: n$_x$="+str(np.round(n.sellmeier_x(lambd1), 4))+"; n$_y$="+str(np.round(n.sellmeier_y(lambd1), 4))
         +"; n$_z$="+str(np.round(n.sellmeier_z(lambd1), 4)), bbox=dict(edgecolor='black', facecolor='white'))
ax1.legend(bbox_to_anchor=(1, 1), loc=2)
ax1.set_title("Phase Matching for Doubling 1055 nm in YSO")
ax1.set_xlabel(r"$\phi$ (degrees)")
ax1.set_ylabel(r"$\theta$ (degrees)")
plt.savefig("Phase Matching Angles in YSO.pdf", bbox_inches="tight")
plt.close()
# plt.show()

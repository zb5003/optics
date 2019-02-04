#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   This program plots the coherence length and second harmonic intensity as the incident angle is varied around the
# ideal phase match angle in YSO.  It also plots the variation in an EFISH as the DC field is varied.  The steps proceed
#  as follows:
#
# 1.) Define physical parameters.
#       a.) Calculate intensity based on power and waist.
# 2.) Initialize YSO class from Birefringence.py.
# 3.) Find the principle indices for the fundamental and second harmonic.
# 4.) Find the polar angle appropriate for Type I phase matching for an azimuthal angle of pi/2 (we will be doing Type I
#     phase matching in our crystal.
# 5.) Initialize an instance of Normal_Surface for the fundamental (fast) and second harmonic (slow).
# 6.) Initiate arrays for plot values.
# 7.) Loop to fill those arrays.
# 8.) Plot the arrays.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

from Birefringence import *
import matplotlib.pyplot as plt
import scipy as sp

# 1.)___________________________________________________________________________________________________________________
lambd1 = 1.055 * 10**(-6)  # fundamental (meters)
lambd3 = 0.527 * 10**(-6)  # second harmonic (meters)
omega3 = pc.c / lambd3  # frequency of second harmonic (Hz)
L = 0.004  # length of the crystal (meters)
n = 1.8  # approximate index of refraction
chi3 = 10**(-22)  # third order susceptibility approximation
E = 150 / 0.01

P = 1  # fundamental power (Watts)
waist = 20 * 10**(-6)  # beam waist (meters)
I1 = P / (sp.pi * waist**2)  # fundamental intensity

# 2.)___________________________________________________________________________________________________________________
yso = YSO()

# 3.)___________________________________________________________________________________________________________________
n1_x = yso.sellmeier_x(lambd1 * 10**(6))
n1_y = yso.sellmeier_y(lambd1 * 10**(6))
n1_z = yso.sellmeier_z(lambd1 * 10**(6))

n3_x = yso.sellmeier_x(lambd3 * 10**(6))
n3_y = yso.sellmeier_y(lambd3 * 10**(6))
n3_z = yso.sellmeier_z(lambd3 * 10**(6))

# 4.)___________________________________________________________________________________________________________________
ysoI = TypeI(n1_x, n1_y, n1_z, n3_x, n3_y, n3_z)
phi = 90 * sp.pi / 180
theta = ysoI.from_equator(phi)

# 5.)___________________________________________________________________________________________________________________
n1_normal = Normal_Surface(n1_x, n1_y, n1_z)
n3_normal = Normal_Surface(n3_x, n3_y, n3_z)

# 6.)___________________________________________________________________________________________________________________
points = 500  # number of plot points

angle_var = 5 * sp.pi / 180  # error in each angle (+/-)
phi_var = sp.linspace(phi - angle_var, phi + angle_var, points)
theta_var = sp.linspace(theta - angle_var, theta + angle_var, points)

Lcoh_phi = sp.zeros(points)  # coherence length as phi is varied
Lcoh_theta = sp.zeros(points)  # coherence length as theta is varied
Lcoh_E = sp.zeros(points)  # coherence length as the DC field is varied

I_phi = sp.zeros(points)  # second harmonic as phi is varied
I_theta = sp.zeros(points)  # second harmonic as theta is varied
I_E = sp.zeros(points)  # second harmonic as the DC field is varied

# 7.)___________________________________________________________________________________________________________________
for i in range(points):
    n1_phi = n1_normal.n2(theta, phi_var[i])
    n3_phi = n3_normal.n1(theta, phi_var[i])

    n1_theta = n1_normal.n2(theta_var[i], phi)
    n3_theta = n3_normal.n1(theta_var[i], phi)

    Lcoh_phi[i] = lambd1 * lambd3 / (2 * ((n3_phi - n1_phi) * lambd1))
    Lcoh_theta[i] = lambd1 * lambd3 / (2 * ((n3_theta - n1_theta) * lambd1))

    I_phi[i] = 8 * (chi3 * E)**2 * omega3**2 * I1**2 * L**2 / (n**3 * pc.epsilon0 * pc.c**2) * sp.sinc(L / Lcoh_phi[i])**2
    I_theta[i] = 8 * (chi3 * E)**2 * omega3**2 * I1**2 * L**2 / (n**3 * pc.epsilon0 * pc.c**2) * \
                 sp.sinc(L / Lcoh_theta[i])**2

# 8.)___________________________________________________________________________________________________________________
fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2)
ax1.plot(phi_var * 180 / sp.pi, I_phi)
ax2.plot(theta_var * 180 / sp.pi, I_theta)
# ax1.set_ylim([0, 15])

plt.show()
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   This program plots the wavelength spectrum of a pulse produced by frequency doubling an ultrashort pulse.
#
# 1.) Define the physical parameters.
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
from Birefringence import *
import matplotlib.pyplot as plt
import time
t1 = time.time()

# 1.)___________________________________________________________________________________________________________________
# All units are SI
L = 0.004  # crystal length
theta = sp.pi / 2  # phase matching polar angle
phi = 21.9 * sp.pi / 180  # phase matching azimuthal angle

I1_amp = 10000 / (sp.pi * 0.0002**2)  # roughly the intensity of the orange laser
I1_center = 1064e-9  # not really the orange lasers center, just for convenience (for now)
I1_width = 15e-9

lambd1_step = 0.2e-9
lambd1_range = 10 * I1_width
n_lambd1 = int(2 * lambd1_range / lambd1_step + 1)
lambd1 = sp.linspace(I1_center - lambd1_range, I1_center + lambd1_range, n_lambd1)

lambd3_step = 1e-9
lambd3_width = 3 * I1_width
lambd3_start = (I1_center - lambd3_width) / 2
lambd3_end = (I1_center + lambd3_width) / 2
n_lambd3 = int(4 * lambd3_width / lambd3_step + 1)
lambd3 = sp.linspace(lambd3_start, lambd3_end, n_lambd3)
I3 = sp.zeros(n_lambd3)  # array to hold sfg
I3_shg = sp.zeros(n_lambd1)  # array to hold shg and sh wavelength

# 2.)___________________________________________________________________________________________________________________
typeI = KTP()
I1_profile = typeI.gaussian_pulse(I1_center, I1_width, I1_amp, lambd1)
for i in range(n_lambd3):
    for j in range(n_lambd1):
        lambd2 = lambd1[j] * lambd3[i] / (lambd1[j] - lambd3[i])
        I1 = typeI.gaussian_pulse(I1_center, I1_width, I1_amp, lambd1[j])
        I2 = typeI.gaussian_pulse(I1_center, I1_width, I1_amp, lambd2)
        I3[i] = I3[i] + typeI.sum_frequenc_I(I1, I2, lambd1[j], lambd3[i], L, theta, phi) * lambd1_step

for k in range(n_lambd1):
    I1 = typeI.gaussian_pulse(I1_center, I1_width, I1_amp, lambd1[k])
    I3_shg[k] = typeI.shg_perfect_type_I(I1, lambd1[k], L, theta, phi)

t2 = time.time()
print("This shit took ", t2-t1, " seconds.")

fig = plt.figure()
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(lambd3, I3)
ax1.set_xlim(432e-9, 632e-9)

ax2 = fig.add_subplot(3, 1, 2)
ax2.plot(lambd1, I1_profile)
ax2.set_xlim(964e-9, 1164e-9)

ax3 = fig.add_subplot(3, 1, 3)
ax3.plot(lambd1 / 2, I3_shg)
ax3.set_xlim(432e-9, 632e-9)

plt.show()

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
from astropy.cosmology import Planck13, z_at_value
import pandas as pd
#--------------------------------------------------------------------------------------------
#variables for these equations:
M_sun = 1.99e33 ## g
BHM2 = (10**6) * M_sun
G = 6.67e-8 ## cm^3 g^-1 s^-2
c = 2.99792458e10 ## cm S^-1
a = 3.086e16 ## cm
r = 3.086e18 #cm this is 1pc

## first, is the original rho equation we used. This comes out to 6x(10^-19).
## using this makes out t_df very small (10^-09)
# rho = ((10.**4.) * M_sun) / (2.93799895e55) #g/cm^3
# print(rho)

## this rho equation was taken from the gxys homework. this comes out to (10^+15)
rho_o = (BHM2 * 1.25) / ( 4 * 3.14159265359 * (r**(1.25)))
print(rho_o)
## making this rho (10^-22). Using this makes t_df ~(5.5*(10^-06)). Still too small.
rho = rho_o * (r**(-7/4))
print(rho)

## the only thing that could concieveablly change is using a different total mass
##of stars versus the seondary BH mass that I did use? The value used in the homework
##was close to this mass, so that's why I figured it was fine, but now I'm nt sure.

# v = 10.**7. #cm/s
v = ((G * BHM2) / r)**(1/2) ## this comes out to 1? why is this wrong?
print(v)

lamduh = 10.
e_const = 2.718281828
t_sal = 1.26144e15 #s
yr = 3.17098e-8 # 1 sec = this many years
t_3body = np.random.uniform((10e6), (3e9), 34.)


dataframe = pd.read_csv('glenna.csv')
r_z = np.array(dataframe['redshift'])
newredshift = np.array(dataframe['redshift'])
M1 = np.array(dataframe['mass1']) * M_sun
M2 = np.array(dataframe['mass2']) * M_sun
M_halo = np.array(dataframe['haloMass'])

#
# u = (c**5.) / (G**2.)
# s = (a**4.) / ((M1 * M2) * (M1 + M2))
# print(s)

## Peters and Matthews GW equation for timescale of a system to merge by GW emission:
## note that we need the equation for a (the assumed separation) to compute this:

##dynamical friction.
##The value of the power of v changes if you use the equation (5) or the solid value we guesses(3)
t_df = (v**5) / (4. * 3.14159265359 * (BHM2) * G**2 * lamduh * rho)

df_years = t_df * 3.154e+7
print(df_years) ## comes out to 176 years. need 90 million more, lol
# t_df = 3e15 # sec

t_gw = (5. / 256.) * ((c**5.) / (G**2.)) * ((a**4.) / ((M1 * M2) * (M1 + M2)))

t_merge = t_df + t_3body + t_gw

## this will be the new secondary BH mass:
M_BH2 = (M2 * (e_const ** ((t_merge * .01) / (t_sal)) - 1)) + M2

## how much growth happens per year
M_edd = 2.26 * M2 * 0.1 # (g / yr)


# print(M1)
# print(M2)
# # print(t_3body)
print(t_df)
# # print(t_gw)
# # print(t_merge)
# print(M_BH2)
# # print(M_edd)


# ## make excel sheet with these values (new mass):
# #
# dataframe_out=pd.DataFrame([t_gw,t_merge,M_BH2],columns=['t_gw','t_merge','M_BH2'])
# print(datafrom_out)
# dataframe_out.to_csv('new_masses.csv')


## new redshifts based on above timescales

# for n in range(len(r_z)):
#
#     ##find age at z= whatever...
#     age = cosmo.age(r_z[n])
#     # age13 = cosmo.age(13)
#
#     t3 = (t_3body[n] + t_df) / (3e16)
#
#     #add 1 billion years
#     newage = age + (t3 * u.Gyr)
#
#
#     ## convert newage (the age plus a billion years) back to redshift:
#     newredshift[n] = z_at_value(Planck13.age, newage)
#
#     # print(r_z)
#     # print (age)
#     # print(newage)
#     # print (newredshift)
#
# BHdata = np.column_stack((M1, M2, M_BH2, r_z, newredshift))
# np.savetxt('/Users/Olivia/Desktop/BH_data.csv', BHdata, fmt='%.8e', delimiter=',', header= 'M1, M2, M_BH2, r_z, newredshift')

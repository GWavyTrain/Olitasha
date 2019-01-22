import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
from astropy.cosmology import Planck13, z_at_value
import pandas as pd
#--------------------------------------------------------------------------------------------

## calculate strain for two SMBH of mass 10^6 each, at a distance of z=.01 and znew=13

## constants:
M_sun = 1.99 * (10.**33.) ## g
# M1 = (10.**6.) * M_sun ## M_sun
# M2 = (10.**6.) * M_sun ## M_sun
# r_znew = 1.667* (10.**29.) ## cm
# r_z = 1.303* (10.**26.) ## cm
G = 6.67 * (10.**-8.) ## cm^3 g^-1 s^-2
c = 2.99792458 * (10.**10.) ## cm S^-1
# f = (10**-3) ## assumed from LISA GW strain/frequency plot
f = np.logspace(np.log10(10e-5), np.log10(1e0), 1000.)


# #--------------------------------------------------------------------------------------------
#
## here is the astropy cosmo module stuff:

dataframe = pd.read_csv('glenna.csv')
r_z = np.array(dataframe['redshift'])
M1 = np.array(dataframe['mass1']) * M_sun
M2 = np.array(dataframe['mass2']) * M_sun

for n in range(len(r_z)):

    ## find age at z= whatever...
    age = cosmo.age(r_z[n])
    # age13 = cosmo.age(13)

    # add 1 billion years
    newage = age + (1 * u.Gyr)

    ## convert newage (the age plus a billion years) back to redshift:
    redshift = z_at_value(Planck13.age, newage)


    # print (age)
    # # print (age13)
    # print (redshift)

    ## these are the calculated values that are now input into the constants section:
    ## @ z=1: age = 5.92222860049 Gyr
    ## add 1 Gyr years to age at z=1 and get: 6.92222860049 Gyr
    ## new redshift with age = 6.9 gyr: 0.770882891415


    #--------------------------------------------------------------------------------------------
    #---------------------------------------------------------------------------------------------

    ## calulate chirp mass:

    M_c = (((M1[n] * M2[n])**(3./5.)) / ((M1[n]+M2[n])**(1./5.)))

    ## calculate pure strain:

    h_pure = (((8. * np.pi**(2./3.)) / (10.**(1./2.))) * (((G**(5./3.)) * (M_c**(5./3.))) / ((c**4.) * (r_z[n])))) * (f**(2./3.))

    h_purenew = (((8. * np.pi**(2./3.)) / (10.**(1./2.))) * (((G**(5./3.)) * (M_c**(5./3.))) / ((c**4.) * (redshift)))) * (f**(2./3.))

    ## Calculate f/f-dot (n):

    n = ((5. / (96.* (np.pi**(8./3.)))) * (c**5. / ((G**(5./3.)) * M_c**(5./3.)))) * f**(-5./3.)

    ## Characteristic Strain:

    hc = h_pure * (n**(1./2.))

    hcnew = h_purenew * (n**(1./2.))

    #--------------------------------------------------------------------------------------------
    # #
    # print (M_c)
    # print (h_pure)
    # print (h_purenew)
    # print (n)
    # print (hc)
    # print (hcnew)

    ## these are the values I got:
    ##  = chirp mass
    ##  = pure strain
    ##  = # n
    ##  = characteristic strain

    ## export values into excel sheet.

    mergerdata = np.column_stack((h_pure, h_purenew, n, hc, hcnew))
    np.savetxt('/Users/Olivia/Desktop/merger_data_'+str(redshift)+'.csv', mergerdata, fmt='%.7e', delimiter=',', header= 'h_pure, h_purenew, n, hc, hcnew')


#--------------------------------------------------------------------------------------------

## plot the characteristic strain against a range of frequencies:
plt.yscale("log")
plt.xscale("log")
#
# plt.plot(f, hc)
# plt.plot(f, hcnew)
# plt.show()


# #variables for these equations:
# #
# a = 3.086*(10**21) #cm
# rho = (10**4) * M_sun #solar masses/pc^3
# v = 10**7 #cm/s
# lamduh = 10
# e = 2.718281828
# t_sal = 40. * (10.**6.)
# t_3body = np.logspace(np.log10(10e9), np.log10(10e6), 1000.)
#
# ## Peters and Matthews GW equation for timescale of a system to merge by GW emission:
# ## note that we need the equation for a (the assumed separation) to compute this:
#
# t_df = v**3 / (4 * np.pi * (10**6) * G * lamduh * rho)
#
# t_gw = (5. / 256) * ((c**5.) / ((G**2.) * (a**4.))) / ((M1 * M2) * (M1 + M2))
# print (t_gw)
#
# t_merge = t_df + t_3body + t_gw
#
# ## this will be the new secondary BH mass:
# M_BH2 = M1 * (e ** ((-t_merge) / (t_sal)))
#
# print(t_df)
# print(t_gw)
# print(t_merge)
# print(M_BH2)

# ## make excel sheet with these values (new mass):
#
# BHdata = np.column_stack((M_BH2))
# np.savetxt('/Users/Olivia/Desktop/BH_data.csv', BHdata, fmt='%.7e', delimiter=',', header= 'M_BH2')
# # --------------------------------------------------------------------------------------------


# ## here we will do time delay for merging:
#
# def get_delay()
#
#     t_delay = 1
#
# return t_delay

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
from astropy.cosmology import Planck13, z_at_value
import pandas as pd

dataframe = pd.read_csv('glenna.csv')
r_z = np.array(dataframe['redshift'])

for n in range(len(r_z)):

    ##find age at z= whatever...
    age = cosmo.age(r_z[n])
    # age13 = cosmo.age(13)

    #add 1 billion years
    newage = age + (1 * u.Gyr)

    ## convert newage (the age plus a billion years) back to redshift:
    newredshift = z_at_value(Planck13.age, newage)


    print (age)
    print (newredshift)

    redshiftdata = np.column_stack((r_z, age, newredshift))
    np.savetxt('/Users/Olivia/Desktop/redshiftdata.csv', redshiftdata, fmt='%.11', delimiter=',', header= 'r_z, age, redshift')

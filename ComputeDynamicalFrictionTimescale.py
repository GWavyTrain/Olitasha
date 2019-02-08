#!/usr/local/bin/python3

import numpy as np
from sys import exit

# Physical constants
G                    = 6.67e-8       # [ cm^3/g / s^2 ]
M_sun                = 1.99e33       # [ g ]
CentimetersPerParsec = 3.086e18      # [ cm ]
YearsPerSecond       = 3.17098e-8    # [ yr ]

# Problem parameters
M_IMBH    = 100.0 * M_sun              # [ g ]
M_SMBH    = 1.0e6 * M_sun              # [ g ]
r         = 1.0 * CentimetersPerParsec # [ cm ]
LogLambda = 10.0                       # [ dimensionless ]

v   = np.sqrt( G * M_SMBH / r ) # [ cm / s ]
rho = ( 1.0e4 * M_sun ) / r**3  # [ g / cm^3 ]

t_df = v**3 / ( 4.0 * np.pi * G**2 * M_IMBH * rho * LogLambda )
print( 'Dynamical friction timescale: {:.2f} Myr'.format( \
        t_df * YearsPerSecond / 1.0e6 ) )

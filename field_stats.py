#!/usr/bin/python python3

import os
import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table

from astropy import units as u
from astropy.coordinates import (SkyCoord, Distance, Galactic, 
                                 EarthLocation, AltAz)
import astropy.coordinates as coord

import matplotlib.pyplot as plt
plt.style.use('presentation.mplstyle')
from astropy.io import fits

slits = Table.read('slits_reduced_z_field.fits')

field = []

grouped = slits.group_by('FIELD')
average_ra = grouped.groups.aggregate(np.mean)
print(average_ra)

for i in range(len(average_ra)):
    c = SkyCoord(ra=average_ra['RA'][i]*u.deg, dec=average_ra['DEC'][i]*u.deg)
    print(c.ra.to_string(u.hourangle, sep=':'))
    print(c.dec.to_string(sep=':'))
                
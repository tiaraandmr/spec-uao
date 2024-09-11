#!/usr/bin/python python3

import os
import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt
plt.style.use('presentation.mplstyle')
from astropy.io import fits

#define path to the .FITS file
home = '/Users/saraswati/Documents/Work/spec-uao/UAO-S118-23A/'
exclude = ["BD", "Feige"]
folder = sorted((f for f in os.listdir(home) if not f.startswith(tuple(exclude))), key=str.lower)

slits_z = Table.read('slits_reduced_z.fits')

field = pd.Series([], dtype='int')

for i in range(len(folder)):

    #list the folder inside home directory
    reduced = home+folder[i]+'/reduced/'
    spectra = sorted(os.listdir(reduced))

    if folder[i] == 'Hviding2023A_0_733':
        field_num = 0
    elif folder[i] == 'Hviding2023A_1_734':
        field_num = 1
    elif folder[i] == 'Hviding2023A_2_735':
        field_num = 2
    elif folder[i] == 'Hviding2023A_4_737':
        field_num = 4
    else:
        continue

    #list the .FITS file for each folder
    prefixes = ["slitA001_1", "slitA015_1", "slitB002_1", "slitB011_1"]
    list_spec = sorted((f for f in os.listdir(reduced+spectra[0]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)

    for spec in list_spec:
        
        # Load the .FITS file
        hdul = fits.open(reduced+spectra[0]+'/obj_abs_1D/'+spec)
        name = hdul[0].header['OBJECT']
        name_clean = name.rsplit('_', 1)[0]

        for j in range(len(slits_z)):
            if str(slits_z['ID'][j]) == str(name_clean):
                field[j] = field_num
            else:
                pass

field = field.sort_index()
slits_z.add_column(field, name='FIELD')
slits_z.write('slits_reduced_z_field.fits', overwrite=True)


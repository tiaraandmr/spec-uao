#!/usr/bin/python python3

import os, glob, shutil, subprocess as sub
from pathlib import Path
import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table, join

import matplotlib.pyplot as plt
from spectres import spectres
from astropy.io import fits

#define path to the .FITS file
home = '/Users/saraswati/Documents/Work/spec-uao/UAO-S118-23A/'
exclude = ["BD", "Feige"]
folder = sorted((f for f in os.listdir(home) if not f.startswith(tuple(exclude))), key=str.lower)

clean_spec = '/Users/saraswati/Documents/Work/spec-uao/clean_spec/'
if os.path.exists(clean_spec):
    shutil.rmtree(clean_spec)
os.makedirs(clean_spec)

for i in range(len(folder)):
    #list the folder inside home directory
    reduced = home+folder[i]+'/reduced/'
    spectra = sorted(os.listdir(reduced))
    print(len(spectra))

    #for folder with 2 days of observation data
    if len(spectra) == 2:
        #list the .FITS file for each observation day
        prefixes = ["slitA001_1", "slitA015_1", "slitB002_1", "slitB011_1"]
        day_1 = sorted((f for f in os.listdir(reduced+spectra[0]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)
        day_2 = sorted((f for f in os.listdir(reduced+spectra[1]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)

        #open .FITS file
        for i in range(len(day_1)):
            #print(reduced+spectra[0]+'/obj_abs_1D/'+day_1[i])

            hdul_1 = fits.open(reduced+spectra[0]+'/obj_abs_1D/'+day_1[i])
            hdul_2 = fits.open(reduced+spectra[1]+'/obj_abs_1D/'+day_2[i])

            #.FITS data
            data_1 = hdul_1[0].data
            data_2 = hdul_2[0].data

            #.FITS header
            header_1 = hdul_1[0].header
            header_2 = hdul_2[0].header

            #object name and pixel properties
            name = header_1['OBJECT']
            name_clean = name.rsplit('_', 1)[0]
            crval1 = header_1['CRVAL1']
            crval2 = header_2['CRVAL1']
            cdelt1 = header_1['CDELT1']

            #flux array
            flux_1 = data_1[0]
            flux_err_1 = data_1[1]

            flux_2 = data_2[0]
            flux_err_2 = data_2[1]

            #mask for bad data
            mask_1 = data_1[3]
            mask_2 = data_2[3]

            #sky lines
            sky_1 = data_1[2]
            sky_2 = data_2[2]

            #build the wavelength array
            wav_1 = crval1 + np.arange(len(flux_1)) * cdelt1
            wav_2 = crval1 + np.arange(len(flux_2)) * cdelt1

            #masking bad data and sky lines
            mask_flux_1 = np.logical_and.reduce([mask_1,sky_1 > 8.0])
            mask_flux_2 = np.logical_and.reduce([mask_2,sky_2 > 8.0])
            flux_1[mask_flux_1] = np.nan
            flux_2[mask_flux_2] = np.nan

            #masking the wavelength with calibration artifact
            lam_1_ma_art = ma.masked_inside(wav_1, 5570, 5585)
            lam_2_ma_art = ma.masked_inside(wav_2, 5570, 5585)

            #create a new wavelength array for both spectra
            new_wav = np.arange(-1, max(len(flux_1),len(flux_2)) + 1) * cdelt1 + (crval1 + crval2) / 2

            #resample both spectra based on the new wavelength array
            flux_1_resample, flux_err_1_resample = spectres(new_wav, lam_1_ma_art, flux_1, spec_errs=flux_err_1, verbose=False)
            flux_2_resample, flux_err_2_resample = spectres(new_wav, lam_2_ma_art, flux_2, spec_errs=flux_err_2, verbose=False)
            new_wav_ma = ma.masked_inside(new_wav, 5570, 5585)
            flux = flux_1_resample/(flux_err_1_resample) + flux_2_resample/(flux_err_2_resample) / (1/flux_err_1_resample**2 + 1/flux_err_2_resample**2)
            flux_err = np.sqrt((flux_err_1_resample**2 + flux_err_2_resample**2)) 
            
            #save spectra into a DataFrame
            spec_new = pd.DataFrame()
            spec_new.insert(0, "WAVELENGTH", new_wav_ma)
            spec_new.insert(1, "FLUX", flux)
            spec_new.insert(2, "FLUX_ERR", flux_err)

            #convert DataFrame into Table
            spec_tab = Table.from_pandas(spec_new)

            #save each spectrum into clean_spec directory
            os.chdir(clean_spec)
            spec_tab.write(name_clean+'.fits')
    
    #for folder with 1 days of observation data
    elif len(spectra) == 1:
        #list the .FITS file for each observation day
        prefixes = ["slitA001_1", "slitA015_1", "slitB002_1", "slitB011_1"]
        day_1 = sorted((f for f in os.listdir(reduced+spectra[0]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)

        #open .FITS file
        for i in range(len(day_1)):
            #print(reduced+spectra[0]+'/obj_abs_1D/'+day_1[i])
            hdul_1 = fits.open(reduced+spectra[0]+'/obj_abs_1D/'+day_1[i])

            #.FITS data
            data_1 = hdul_1[0].data

            #.FITS header
            header_1 = hdul_1[0].header

            #object name and pixel properties
            name = header_1['OBJECT']
            name_clean = name.rsplit('_', 1)[0]
            crval1 = header_1['CRVAL1']
            cdelt1 = header_1['CDELT1']

            #flux array
            flux_1 = data_1[0]
            flux_err_1 = data_1[1]

            #mask for bad data
            mask_1 = data_1[3]

            #sky lines
            sky_1 = data_1[2]

            #build the wavelength array
            wav_1 = crval1 + np.arange(len(flux_1)) * cdelt1

            #masking bad data and sky lines
            mask_flux_1 = np.logical_and.reduce([mask_1,sky_1 > 8.0])
            flux_1[mask_flux_1] = np.nan

            #masking the wavelength with calibration artifact
            lam_1_ma_art = ma.masked_inside(wav_1, 5570, 5585)

            #final flux and flux_err
            flux = flux_1
            flux_err = flux_err_1

            #save spectra into a DataFrame
            spec_new = pd.DataFrame()
            spec_new.insert(0, "WAVELENGTH", lam_1_ma_art)
            spec_new.insert(1, "FLUX", flux)
            spec_new.insert(2, "FLUX_ERR", flux_err)

            #convert DataFrame into Table
            spec_tab = Table.from_pandas(spec_new)

            #save each spectrum into clean_spec directory
            os.chdir(clean_spec)
            spec_tab.write(name_clean+'.fits')
    
    #for folder with no observation data
    elif len(spectra) == 0:
        print('No Observation Data')

#list filename inside the clean_spec directory and assign an int type
spectra = glob.glob(clean_spec+'/*')
IDs = np.array([Path(s).stem for s in spectra]).astype(int)

#join the new object list with the catalog list, then create a new table                 
table1 = Table.read('/Users/saraswati/Documents/Work/spec-uao/slits.fits')
table2 = Table([IDs],names=("ID",))
new_table = join(table1,table2) 
new_table.write('/Users/saraswati/Documents/Work/spec-uao/slits_reduced.fits', overwrite=True)
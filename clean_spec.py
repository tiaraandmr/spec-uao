#!/usr/bin/python python3

import os, subprocess as sub
from pathlib import Path
import numpy as np
import numpy.ma as ma
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

from spectres import spectres
from specutils import Spectrum1D
from astropy.io import fits

#define path to the .FITS file
home = '/Users/saraswati/Documents/Work/spec-uao/UAO-S118-23A/'
exclude = ["BD", "Feige"]
folder = sorted((f for f in os.listdir(home) if not f.startswith(tuple(exclude))), key=str.lower)


for i in range(len(folder)):
    #list the folder inside home directory
    reduced = home+folder[i]+'/reduced/'
    spectra = sorted(os.listdir(reduced))
    print(len(spectra))

    if len(spectra) == 2:
        #list the .FITS file for each observation day
        prefixes = ["slitA001_1"]
        day_1 = sorted((f for f in os.listdir(reduced+spectra[0]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)
        day_2 = sorted((f for f in os.listdir(reduced+spectra[1]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)

        #open .FITS file
        for i in range(len(day_1)):
            print(reduced+spectra[0]+'/obj_abs_1D/'+day_1[i])
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

            #automatically set the x-limit for every spectra
            good_1 = np.invert(np.isnan(flux_1))
            good_2 = np.invert(np.isnan(flux_2))

            #create a new wavelength array for both spectra
            new_wav = np.arange(-1, max(len(flux_1),len(flux_2)) + 1) * cdelt1 + (crval1 + crval2) / 2

            #resample both spectra based on the new wavelength array
            flux_1_resample, flux_err_1_resample = spectres(new_wav, lam_1_ma_art, flux_1, spec_errs=flux_err_1, verbose=False)
            flux_2_resample, flux_err_2_resample = spectres(new_wav, lam_2_ma_art, flux_2, spec_errs=flux_err_2, verbose=False)
            total_flux = flux_1_resample/(flux_err_1_resample) + flux_2_resample/(flux_err_2_resample) / (1/flux_err_1_resample**2 + 1/flux_err_2_resample**2)
            total_flux_err = np.sqrt((flux_err_1_resample**2 + flux_err_2_resample**2)) 
    
    elif len(spectra) == 1:
        #list the .FITS file for each observation day
        prefixes = ["slitA001_1"]
        day_1 = sorted((f for f in os.listdir(reduced+spectra[0]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)

        #open .FITS file
        for i in range(len(day_1)):
            print(reduced+spectra[0]+'/obj_abs_1D/'+day_1[i])
            hdul_1 = fits.open(reduced+spectra[0]+'/obj_abs_1D/'+day_1[i])

            #.FITS data
            data_1 = hdul_1[0].data

            #.FITS header
            header_1 = hdul_1[0].header

            #object name and pixel properties
            name = header_1['OBJECT']
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

            #automatically set the x-limit for every spectra
            good_1 = np.invert(np.isnan(flux_1))

            total_flux = flux_1
            total_err = flux_err_1
    
    elif len(spectra) == 0:
        print('No Observation Data')
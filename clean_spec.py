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
from astropy.wcs import WCS

#delete .DS_Store file in every directory inside spec-uao if there's any
sub.run(["find", ".", "-name ", ".DS_Store", "-delete"])

#define path to the .FITS file
home = '/Users/saraswati/Documents/Work/spec-uao/UAO-S118-23A/'
folder = home+'Hviding2023A_0_733/reduced/'

#list the folder inside home directory
spectra = os.listdir(folder)

#list the .FITS file for each observation day
for i in range(len(spectra)):
    day_1 = sorted(os.listdir(folder+spectra[0]+'/obj_abs_1D/'))
    day_2 = sorted(os.listdir(folder+spectra[1]+'/obj_abs_1D/'))

    #open .FITS file
    for i in range(len(day_1)):
        hdul_1 = fits.open(folder+spectra[0]+'/obj_abs_1D/'+day_1[i])
        hdul_2 = fits.open(folder+spectra[0]+'/obj_abs_1D/'+day_2[i])

        #.FITS data
        data_1 = hdul_1[0].data
        data_2 = hdul_2[0].data

        #.FITS header
        header_1 = hdul_1[0].header
        header_2 = hdul_2[0].header

        #object name and pixel properties
        name = header_1['OBJECT']
        crval1 = header_1['CRVAL1']
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
        w_1 = WCS(header_1, naxis=1, relax=False, fix=False)
        lam_1 = (w_1.wcs_pix2world(np.arange(len(flux_1)), 0)[0])*10**10

        w_2 = WCS(header_2, naxis=1, relax=False, fix=False)
        lam_2 = (w_2.wcs_pix2world(np.arange(len(flux_2)), 0)[0])*10**10

        #masking the bad data
        flux_1_ma = ma.masked_array(flux_1, mask_1)
        flux_2_ma = ma.masked_array(flux_2, mask_2)

        #masking the sky lines
        flux_1_ma_sky = ma.masked_greater(flux_1_ma, 8.0)
        flux_2_ma_sky = ma.masked_greater(flux_2_ma, 8.0)

        #masking the wavelength with calibration artifact
        lam_1_ma_art = ma.masked_inside(lam_1, 5570, 5585)
        lam_2_ma_art = ma.masked_inside(lam_2, 5570, 5585)

        #automatically set the x-limit for every spectra
        good_1 = np.invert(np.isnan(flux_1))
        good_2 = np.invert(np.isnan(flux_2))
        
        #resample second day spectra
        regrid = np.arange(crval1, lam_1[good_1].max(), cdelt1) + 2.5
        spec_resample, spec_errs_resample = spectres(regrid, lam_2_ma_art, flux_2_ma_sky, spec_errs=flux_err_2)
        lam_2_resample = ma.masked_inside(regrid, 5570, 5585)
        flux_2_resample = spec_resample
        flux_err_2_resample = spec_errs_resample

        #flux_total = (flux_1_ma_sky + flux_2_resample) / 2
        #err_total = np.sqrt(flux_err_1**2 + flux_err_2_resample**2)
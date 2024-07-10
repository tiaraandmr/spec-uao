#!/usr/bin/python python3

import os, shutil, subprocess as sub
from pathlib import Path
import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table
from spectres import spectres


import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

from astropy.io import fits

day_1 = '/Users/saraswati/Documents/Work/spec-uao/UAO-S118-23A/Hviding2023A_1_734/reduced/2023.0412/obj_abs_1D/slitA002_41231978099401625t95686.fits'
day_2 = '/Users/saraswati/Documents/Work/spec-uao/UAO-S118-23A/Hviding2023A_1_734/reduced/2023.0415/obj_abs_1D/slitA002_41231978099401625t95686.fits'

hdul_1 = fits.open(day_1)
hdul_2 = fits.open(day_2)

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
mask_flux_1 = np.logical_and(mask_1,sky_1)
mask_flux_2 = np.logical_and(mask_2,sky_2)
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

#automatically set the x-limit for every spectra
good = np.invert(np.isnan(flux_1))
good_2 = np.invert(np.isnan(flux_2))

#define the figure and font size
plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(2, figsize=(17,5), sharex=True, gridspec_kw={'height_ratios': [3,1]})
fig.subplots_adjust(hspace=0, wspace=0)

#plot the spectrum and set the x-limit
axes[0].plot(new_wav_ma, flux_1_resample, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Night 1 (N1)')
axes[0].plot(new_wav_ma, flux_2_resample, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label = 'Night 2 (N2)', alpha=0.5)
axes[0].plot(wav_2, sky_2, color='black', linewidth=2.0, drawstyle='steps-mid')
axes[0].set_xlim(wav_1[good].min(), wav_1[good].max())

#plot the difference between 2 spectrum
delta = flux_1_resample - flux_2_resample 
axes[1].plot(new_wav_ma, delta, color='black', linewidth=2.0, drawstyle='steps-mid')

#define x and y label and plot title
axes[1].set_xlabel(r'Observed Wavelength [$ \rm \AA$]')
axes[0].set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]')
axes[1].set_ylabel(r'N1--N2')
title = Path(day_1).stem
fig.suptitle(title)
axes[0].legend()

#save the plot as .pdf
fig.savefig(title+'.pdf', bbox_inches="tight")
plt.close(fig)
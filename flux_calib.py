#!/usr/bin/python python3

import os, glob, shutil, subprocess as sub
from pathlib import Path
import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt
from spectres import spectres
from astropy.io import fits
import matplotlib.ticker as ticker
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

home = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/'

feige_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/ffeige34.dat'
feige = np.genfromtxt(feige_path)

feige_path_obs = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/clean_calib/Feige34_B_clean.fits'
feige_obs = Table.read(feige_path_obs)

feige_path_orig = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/Feige34_B.fits'
hdul = fits.open(feige_path_orig)

#resample the template spectrum to the wavelength range of the observed spectrum
cdelt1 = hdul[0].header['CDELT1']
crval1 = hdul[0].header['CRVAL1']

flux_1 = hdul[0].data

#create a new wavelength array for both spectra
new_wav = np.arange(-1, len(feige_obs['FLUX']) + 1) * cdelt1 + crval1

#resample both spectra based on the new wavelength array
flux_feige_resample, flux_err_feige_resample = spectres(new_wav, feige[:,0], feige[:,1]*10, spec_errs=feige[:,2], verbose=False)
flux_obs_resample, flux_err_obs_resample = spectres(new_wav, feige_obs['WAVELENGTH'], feige_obs['FLUX'], spec_errs=feige_obs['FLUX_ERR'], verbose=False)

#calculate the flux sampling
sampling = flux_feige_resample / flux_obs_resample 

#automatically set the x-limit for every spectra
good = np.invert(np.isnan(sampling))

#define the figure and font size
plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(3, figsize=(17,5), sharex=True, gridspec_kw={'height_ratios': [3,3,1.5]})
fig.subplots_adjust(hspace=0, wspace=0)

#plot the spectrum and set the x-limit
axes[0].plot(new_wav, flux_feige_resample, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Template')
#axes[0].plot(feige[:,0], feige[:,1], color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Template')
axes[0].legend()
axes[1].plot(new_wav, flux_obs_resample, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Observed', alpha=0.7)
#axes[1].plot(feige_obs['WAVELENGTH'], feige_obs['FLUX'], color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Observed')
axes[1].legend()
axes[2].plot(new_wav, sampling, color='black', linewidth=2.0, drawstyle='steps-mid')
#axes[2].set_xlim(6500,9000)
axes[2].set_xlim(new_wav[good].min(), new_wav[good].max())
                
#define x and y label and plot title
axes[2].set_xlabel(r'Observed Wavelength [$ \rm \AA$]')
axes[0].set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]')
axes[0].yaxis.set_label_coords(-0.06,-0.03)
axes[2].set_ylabel(r'F($\lambda$)', labelpad=25)
fig.suptitle('Feige 34')

#save the plot as .pdf
fig.savefig(home+'Feige34.pdf', bbox_inches="tight")
plt.close(fig)
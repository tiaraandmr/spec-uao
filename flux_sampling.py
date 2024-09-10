#!/usr/bin/python python3

# Import packages
import os, shutil
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from astropy.table import Table
from spectres import spectres
import pandas as pd

# Set the font style to LaTeX
rc('font', **{'family': 'serif', 'serif': ['lmodern']})
rc('text', usetex=True)

# Define home directory
home = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/'

# Load the template
calib_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/calib_dat/fbd33d2642.dat'
calib = np.genfromtxt(calib_path)
calib[:, 1:2] *= 10  # Convert to same units as observed spectrum

# Load the observed spectrum
calib_path_obs = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/clean_calib/BD+33d2642_B_20230415_clean.fits'
calib_obs = Table.read(calib_path_obs)

# Create a new wavelength array for both spectra
wav = calib_obs['WAVE']

# Resample both spectra based on the new wavelength array
flux_calib_resample, flux_err_calib_resample = spectres(wav, calib[:, 0], calib[:, 1], spec_errs=calib[:, 2], verbose=False)

# Calculate the flux sampling
sampling = flux_calib_resample / calib_obs['FLUX']
sampling_err = sampling * np.sqrt((flux_err_calib_resample/flux_calib_resample)**2 + (calib_obs['ERROR']/calib_obs['FLUX'])**2)

# Automatically set the x-limit for every spectra
good = np.invert(np.isnan(sampling))

# Save spectra into a DataFrame and convert DataFrame into Table
mask = np.logical_or(np.isnan(sampling), np.isnan(sampling_err))
spec_tab = Table([wav, sampling, sampling_err, mask], names=('WAVE', 'FLUX', 'ERROR', 'MASK'))

# Save each spectrum into clean_spec directory
spec_tab.write(home+'BD+33d2642_20230415_sampling.fits')

# Define the figure and font size
plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(2, figsize=(17, 5), sharex=True, gridspec_kw={'height_ratios': [6, 1.5]})
fig.subplots_adjust(hspace=0, wspace=0)

# Plot the spectrum and set the x-limit
axes[0].plot(wav, flux_calib_resample, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Template', alpha=0.7)
axes[0].plot(wav, calib_obs['FLUX'], color='steelblue', linewidth=2.0, drawstyle='steps-mid', label='Observed', alpha=0.7)
axes[0].legend()
axes[1].plot(wav, sampling, color='black', linewidth=2.0, drawstyle='steps-mid')
axes[1].set_xlim(wav[good].min(), wav[good].max())

# Define x and y label and plot title
axes[1].set_xlabel(r'Observed Wavelength [$ \rm \AA$]')
axes[0].set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]')
axes[1].set_ylabel(r'F($\lambda$)', labelpad=17)
fig.suptitle('BD+33d2642_20230415')

# Save the plot as .pdf
fig.savefig(home+'BD+33d2642_20230415.pdf', bbox_inches='tight')
plt.close(fig)
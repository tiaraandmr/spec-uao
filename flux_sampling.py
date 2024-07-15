#!/usr/bin/python python3

# Import packages
import os
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
calib_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/ffeige34.dat'
calib = np.genfromtxt(calib_path)
calib[:, 1:2] *= 10  # Convert to same units as observed spectrum

# Load the observed spectrum
calib_path_obs = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/clean_calib/Feige34_B_clean.fits'
calib_obs = Table.read(calib_path_obs)

# create a new wavelength array for both spectra
new_wav = calib_obs['WAVELENGTH']

# resample both spectra based on the new wavelength array
flux_calib_resample, flux_err_calib_resample = spectres(new_wav, calib[:, 0], calib[:, 1], spec_errs=calib[:, 2], verbose=False)

# calculate the flux sampling
sampling = flux_calib_resample / calib_obs['FLUX']

# automatically set the x-limit for every spectra
good = np.invert(np.isnan(sampling))

#save spectra into a DataFrame
spec_new = pd.DataFrame()
spec_new.insert(0, "WAVELENGTH", new_wav)
spec_new.insert(1, "FLUX", sampling)

#convert DataFrame into Table
spec_tab = Table.from_pandas(spec_new)

#save each spectrum into clean_spec directory
#os.chdir(home)
spec_tab.write(home+'Feige34_sampling.fits')

# define the figure and font size
plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(2, figsize=(17, 5), sharex=True, gridspec_kw={'height_ratios': [6, 1.5]})
fig.subplots_adjust(hspace=0, wspace=0)

# plot the spectrum and set the x-limit
axes[0].plot(new_wav, flux_calib_resample, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Template', alpha=0.7)
axes[0].plot(new_wav, calib_obs['FLUX'], color='steelblue', linewidth=2.0, drawstyle='steps-mid', label='Observed', alpha=0.7)
axes[0].legend()
axes[1].plot(new_wav, sampling, color='black', linewidth=2.0, drawstyle='steps-mid')
axes[1].set_xlim(new_wav[good].min(), new_wav[good].max())

# define x and y label and plot title
axes[1].set_xlabel(r'Observed Wavelength [$ \rm \AA$]')
axes[0].set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]')
axes[1].set_ylabel(r'F($\lambda$)', labelpad=17)
fig.suptitle('Feige34')

# save the plot as .pdf
fig.savefig(home+'Feige34.pdf', bbox_inches='tight')
plt.close(fig)
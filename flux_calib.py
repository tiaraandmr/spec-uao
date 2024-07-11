#!/usr/bin/python python3

import os, glob, shutil, subprocess as sub
from pathlib import Path
import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt
from spectres import spectres
from astropy.io import fits, ascii
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

home = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/'

feige_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/ffeige34.dat'
feige = np.genfromtxt(feige_path)

feige_path_obs = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/clean_calib/Feige34_B_clean.fits'
feige_obs = Table.read(feige_path_obs)

#automatically set the x-limit for every spectra
good = np.invert(np.isnan(feige[:,1]))
good_obs = np.invert(np.isnan(feige_obs['FLUX']))

#resample the template spectrum to the wavelength range of the observed spectrum

#define the figure and font size
plt.rcParams.update({'font.size': 14})
fig, axes = plt.subplots(3, figsize=(17,5), sharex=True, gridspec_kw={'height_ratios': [2,2,1]})
fig.subplots_adjust(hspace=0, wspace=0)

#plot the spectrum and set the x-limit
axes[0].plot(feige[:,0], feige[:,1], color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Template')
axes[0].legend()
axes[1].plot(feige_obs['WAVELENGTH'], feige_obs['FLUX'], color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Observed')
axes[1].legend()
axes[2].plot()
axes[2].set_xlim(feige[:,0][good].min(), feige[:,0][good].max())
                
#define x and y label and plot title
axes[2].set_xlabel(r'Observed Wavelength [$ \rm \AA$]')
axes[1].set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]')
fig.suptitle('Feige 34')

#save the plot as .pdf
fig.savefig(home+'Feige34.pdf', bbox_inches="tight")
plt.close(fig)
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
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

#define path to the .FITS file
home = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/'
exclude = ["plot_calib", "clean_calib", "calib_dat"]
calib = sorted((f for f in os.listdir(home) if not f.startswith(tuple(exclude))), key=str.lower)

plot_calib = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/plot_calib/'
if os.path.exists(plot_calib):
    shutil.rmtree(plot_calib)
os.makedirs(plot_calib)

clean_calib = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/clean_calib/'
if os.path.exists(clean_calib):
    shutil.rmtree(clean_calib)
os.makedirs(clean_calib)

for i in range(len(calib)):
    #open .FITS file
    hdul = fits.open(home+calib[i])

    #.FITS data
    data = hdul[0].data

    #.FITS header
    header = hdul[0].header

    #object name and pixel properties
    name = header['OBJECT']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']

    #flux array
    flux = data[0]
    flux_err = data[1]

    #mask for bad data
    mask = data[3]

    #sky lines
    sky = data[2]

    #build the wavelength array
    wav = crval1 + np.arange(len(flux)) * cdelt1

    #masking bad data and sky lines
    mask_flux = np.logical_and(mask, sky)
    flux[mask_flux] = np.nan
    flux_err[mask_flux] = np.nan

    #masking the wavelength with calibration artifact
    lam_ma_art = ma.masked_inside(wav, 5570, 5585)
    lam_ma_art_2 = ma.masked_inside(wav, 6470, 6490)

    #save spectra into a DataFrame
    spec_new = pd.DataFrame()
    spec_new.insert(0, "WAVELENGTH", lam_ma_art_2)
    spec_new.insert(1, "FLUX", flux)
    spec_new.insert(2, "ERR", flux_err)

    #convert DataFrame into Table
    spec_tab = Table.from_pandas(spec_new)

    #save each spectrum into clean_spec directory
    title = Path(calib[i]).stem
    os.chdir(clean_calib)
    spec_tab.write(title+'_clean.fits')

    #automatically set the x-limit for every spectra
    good = np.invert(np.isnan(flux))

    #define the figure and font size
    plt.figure(figsize=(17,5))
    plt.rcParams.update({'font.size': 17})

    #plot the spectrum and set the x-limit
    plt.plot(lam_ma_art_2, flux, color='firebrick', linewidth=2.0, drawstyle='steps-mid')
    plt.xlim(wav[good].min(), wav[good].max())
                
    #define x and y label and plot title
    plt.xlabel(r'Observed Wavelength [$ \rm \AA$]')
    plt.ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]')
    plt.title(title)

    #save the plot as .pdf
    plt.savefig(plot_calib+title+'.pdf', bbox_inches="tight")
    plt.close()
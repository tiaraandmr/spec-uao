#!/usr/bin/python python3

# Import packages
import os, shutil
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from astropy.table import Table
from spectres import spectres
from astropy.convolution import convolve

from astropy.io import fits
import pandas as pd

# Set the font style to LaTeX
rc('font', **{'family': 'serif', 'serif': ['lmodern']})
rc('text', usetex=True)

#define path to the .FITS file
home = '/Users/saraswati/Documents/Work/spec-uao/UAO-S118-23A/'
exclude = ["BD", "Feige"]
folder = sorted((f for f in os.listdir(home) if not f.startswith(tuple(exclude))), key=str.lower)

calibrated_spec = '/Users/saraswati/Documents/Work/spec-uao/calibrated_spec/'
if os.path.exists(calibrated_spec):
    shutil.rmtree(calibrated_spec)
os.makedirs(calibrated_spec)

for f in folder:

    #list the folder inside home directory
    reduced = home+f+'/reduced/'
    spectra = sorted(os.listdir(reduced))

    if f == 'Hviding2023A_0_733':
        sampling_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/fits_sampling/Feige34_sampling.fits'
    elif f == 'Hviding2023A_1_734':
        sampling_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/fits_sampling/BD+33d2642_20230415_sampling.fits'
    elif f == 'Hviding2023A_2_735':
        sampling_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/fits_sampling/BD+33d2642_20230708_sampling.fits'
    elif f == 'Hviding2023A_4_737':
        sampling_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/fits_sampling/BD+33d2642_20230708_sampling.fits'
    else:
        continue

    #list the .FITS file for each folder
    prefixes = ["slitA001_1", "slitA015_1", "slitB002_1", "slitB011_1"]
    list_spec = sorted((f for f in os.listdir(reduced+spectra[0]+'/obj_abs_1D/') if not f.startswith(tuple(prefixes))), key=str.lower)

    for spec in list_spec:

        # Load the flux sampling template
        sampling = Table.read(sampling_path)

        # Load the .FITS file
        hdul = fits.open(reduced+spectra[0]+'/obj_abs_1D/'+spec)
        name = hdul[0].header['OBJECT']
        name_clean = name.rsplit('_', 1)[0]

        # Load the observed spectrum
        spec_obs_path = '/Users/saraswati/Documents/Work/spec-uao/clean_spec/'+name_clean+'.fits'
        spec_obs = Table.read(spec_obs_path)

        # create a new wavelength array for both spectra
        new_wav = spec_obs['WAVE']
        flux = spec_obs['FLUX']
        err = spec_obs['ERROR']
        mask = spec_obs['MASK']

        # resample the 'sampling' spectra based on the new wavelength array
        flux_sampling_resample, flux_err_sampling_resample = spectres(new_wav, sampling['WAVELENGTH'], sampling['FLUX'], spec_errs=sampling['ERR'], verbose=False)

        # Apply the flux sampling to the observed spectra
        flux_calibrated = spec_obs['FLUX'] * flux_sampling_resample
        flux_err_calibrated = flux_calibrated * np.sqrt((flux_err_sampling_resample/flux_sampling_resample)**2 + (spec_obs['ERROR']/spec_obs['FLUX'])**2)

        #smooth the spectra
        N = 5
        # spec_convolve = convolve(flux_calibrated, np.ones(N), boundary='extend', normalize_kernel=True)
        # spec_err_convolve = convolve(flux_err_calibrated, np.ones(N), boundary='extend', normalize_kernel=True)

        #save spectra into a DataFrame
        spec_new = pd.DataFrame()
        spec_new.insert(0, "WAVE", new_wav)
        spec_new.insert(1, "FLUX", flux_calibrated)
        spec_new.insert(2, "ERROR", flux_err_calibrated)
        spec_new.insert(3, "MASK", mask)

        #convert DataFrame into Table
        spec_tab = Table.from_pandas(spec_new)

        #save each spectrum into clean_spec directory
        os.chdir(calibrated_spec)
        spec_tab.write(name_clean+'.fits')
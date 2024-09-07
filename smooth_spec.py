#!/usr/bin/python3

import os, glob, shutil, subprocess as sub
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.convolution import convolve

from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
from astropy.io import fits

# Define source folder
clean_spec = '/home/tiara/Documents/Work/spec-uao/clean_spec/'

# Create new folder for the results
smooth_spec = '/home/tiara/Documents/Work/spec-uao/smooth_spec/'
if os.path.exists(smooth_spec):
    shutil.rmtree(smooth_spec)
os.makedirs(smooth_spec)

# Define line species
species = np.array([['[SII]', 6716.435709204644], ['[SII]', 6730.811841618589], ['[NII]', 6583.451474377235], ['[NII]', 6548.050991532048], 
                    ['[OIII]', 5006.843330296288], ['[OIII]', 4958.91106815638], ['[OIII]', 4363.2096874449635], ['[NeIII]', 3868.7632032493057], 
                    ['[NeV]', 3425.8676240424247], ['[NeV]', 3345.828076914398], ['[OI]', 6300.303551730963], ['[OI]', 6363.776506749638], 
                    ['[OII]', 3727.4199178007175], ['HI', 6562.797027356974], ['HI', 4861.321979760415], ['HI', 4340.459677187083], 
                    ['CIV', 1548.8576020874154], ['HeII', 1639.7908206670402], ['CIII]', 1908.1020362243723], ['MgII', 2798.2921038031536]])

# Define object ID and preliminary z
slits = Table.read('/home/tiara/Documents/Work/spec-uao/slits_reduced.fits')
object_ID = slits['ID'].astype('str')
z = slits['BEST_Z']

for i in range(len(object_ID)):
    # Read the spectrum
    spec_data = Table.read(clean_spec+object_ID[i]+'.fits')
    wav = spec_data['WAVE']
    flux = spec_data['FLUX']
    
    # Smooth the spectra with N=5
    N = 5
    spec_convolve = convolve(flux, np.ones(N), boundary='extend', normalize_kernel=True)
    good = np.invert(np.isnan(flux))

    # Define the figure and font size
    plt.rcParams.update({'font.size': 17})
    fig, ax = plt.subplots(figsize=(17,5))

    # Plot the spectrum and set the x-limit
    ax.plot(wav, spec_convolve, color='firebrick', linewidth=2.0, drawstyle='steps-mid')
    
    # Define the approximate redshift (z)
    best_z = z[i]

    #check if lines are in the range of the observed spectrum
    for j in range(len(species)):
        data = species[j][1].astype('float64')
        x = data*(1+best_z)

        if ((x < wav[-1]) and (x > wav[0])):
            ax.vlines(x, ymin=flux.min(), ymax=flux.max(), color = 'black', ls = '--') 
            ax.text(x, 1.25*flux.max(), species[j][0], rotation=90, fontsize=14, ha='center', va='center')

    #setting the x- and y-axis limit
    ax.set_xlim(wav[good].min(), wav[good].max())
    ax.set_ylim(flux.min(), 1.5*flux.max())

    #plot a rest frame wavelength based on the best_z
    ax2 = ax.secondary_xaxis('top',functions=(lambda lam: lam/(1+best_z), lambda lam: lam*(1+best_z))) 
    ax2.set_xlabel(r'Rest Wavelength [$ \rm \AA$]', labelpad=10)
    
    #define x and y label and plot title
    ax.set_xlabel(r'Observed Wavelength [$ \rm \AA$]', labelpad=10)
    ax.set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]', labelpad=10)
    title = object_ID[i]
    ax.set_title(title+r', $z = {:.3f}$'.format(best_z), pad=20)

    #save the plot as .pdf
    plt.savefig(smooth_spec+title+'_smooth_n5.pdf', bbox_inches="tight")
    plt.close()
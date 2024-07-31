#!/usr/bin/python python3

import os, shutil, subprocess as sub
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

from astropy.io import fits

#define directory and list the file inside clean_spec directory
clean_spec = '/Users/saraswati/Documents/Work/spec-uao/clean_spec/'
list_spec = sorted(os.listdir(clean_spec))

plot_spec = '/Users/saraswati/Documents/Work/spec-uao/plot_spec/'
if os.path.exists(plot_spec):
    shutil.rmtree(plot_spec)
os.makedirs(plot_spec)

for i in range(len(list_spec)):
    #read the spectrum
    spec_data = Table.read(clean_spec+list_spec[i])

    #automatically set the x-limit for every spectra
    good = np.invert(np.isnan(spec_data['FLUX']))

    #define the figure and font size
    plt.figure(figsize=(17,5))
    plt.rcParams.update({'font.size': 17})

    #plot the spectrum and set the x-limit
    plt.plot(spec_data['WAVE'], spec_data['FLUX'], color='firebrick', linewidth=2.0, drawstyle='steps-mid')
    plt.xlim(spec_data['WAVE'][good].min(),spec_data['WAVE'][good].max())
    
    #define x and y label and plot title
    plt.xlabel(r'Observed Wavelength [$ \rm \AA$]')
    plt.ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]')
    title = Path(list_spec[i]).stem
    plt.title(title)

    #save the plot as .pdf
    plt.savefig(plot_spec+title+'.pdf', dpi=1000, bbox_inches="tight")
    plt.close()
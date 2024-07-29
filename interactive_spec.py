#!/usr/bin/env python3

# Import packages
import os
import numpy as np
from astropy.table import Table
from astropy.convolution import convolve

from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
from astropy.io import fits
from spectres import spectres

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--id', default='41231836365464984')
args = parser.parse_args()
spec_id = args.id

# Load the observed spectrum
spec_obs_path = f'/Users/saraswati/Documents/Work/spec-uao/clean_spec/{spec_id}.fits'
spec_obs = Table.read(spec_obs_path)

# Create a new wavelength array for both spectra
new_wav = spec_obs['WAVELENGTH']
flux_calibrated = spec_obs['FLUX']

# Automatically set the x-limit for every spectrum
good = np.invert(np.isnan(flux_calibrated))

# Smooth the spectra
N = 5
spec_convolve = convolve(flux_calibrated, np.ones(N), boundary='extend', normalize_kernel=True)
spec_convolve_orig = convolve(spec_obs['FLUX'], np.ones(N), boundary='extend', normalize_kernel=True)

# Define line species
species = np.array([['[SII]', 6716.435709204644], ['[SII]', 6730.811841618589], ['[NII]', 6583.451474377235], ['[NII]', 6548.050991532048], 
                    ['[OIII]', 5006.843330296288], ['[OIII]', 4958.91106815638], ['[OIII]', 4363.2096874449635], ['[NeIII]', 3868.7632032493057], 
                    ['[NeV]', 3425.8676240424247], ['[NeV]', 3345.828076914398], ['[OI]', 6300.303551730963], ['[OI]', 6363.776506749638], 
                    ['[OII]', 3727.4199178007175], ['HI', 6562.797027356974], ['HI', 4861.321979760415], ['HI', 4340.459677187083], 
                    ['CIV', 1548.8576020874154], ['HeII', 1639.7908206670402], ['CIII]', 1908.1020362243723], ['MgII', 2798.2921038031536]])

# define the figure and font size
plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(17, 5))
fig.subplots_adjust(hspace=0, wspace=0, bottom=0.2)

ax.plot(new_wav, spec_convolve, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Calibrated \& Smoothed (N=5)', alpha=0.7)

#define x and y label and plot title
ax.set_xlabel(r'Observed Wavelength [$ \rm \AA$]', labelpad=10)
ax.set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]', labelpad=10)

#setting the x- and y-axis limit
ax.set_xlim(new_wav[good].min(), new_wav[good].max())
ax.set_ylim(np.nanmin(spec_convolve), 1.5*np.nanmax(spec_convolve))

# Plot lines
lines = [ax.vlines(line, np.nanmin(spec_convolve), np.nanmax(spec_convolve), color = 'black', ls = '--')
         for line in species[:, 1].astype('float64')]
labels = [ax.text(float(line),1.25*np.nanmax(spec_convolve),label,ha='center',va='center',rotation=90) for label,line in species]

# Update function
def submit(redshift):
    global z
    z = float(redshift)

    # Update line positions
    for line_obj, line in zip(lines, species[:, 1].astype('float64')):
        line_obj.set_segments(
            [[(line * (1 + z), np.nanmin(spec_convolve)), (line * (1 + z), 1.1*np.nanmax(spec_convolve))]]
        )
    for label, line in zip(labels, species[:, 1].astype('float64')):
        label.set_x(line * (1+z))
    
    ax.relim()
    ax.autoscale_view()
    plt.draw()

axbox = fig.add_axes([0.1, 0.05, 0.8, 0.075])
text_box = TextBox(axbox, 'Redshift', textalignment='center')
text_box.on_submit(submit)
text_box.set_val('0.8')  # Trigger `submit` with the initial string.

plt.show()

#define path to redshift.txt
redshift = '/Users/saraswati/Documents/Work/spec-uao/redshift.txt/'

#write results to redshift.txt
if os.path.exists('redshift.txt'):
    id = np.genfromtxt('redshift.txt', usecols=0, dtype='unicode')
    z_data = np.genfromtxt('redshift.txt', usecols=1, dtype='unicode')

    if spec_id not in id:
        with open('redshift.txt', 'a') as out:
            out.write('\n'+ spec_id + '\t' + str(z))
            print(str(z))
    
    else:
        for i in range(len(z_data)):
            if spec_id == id[i]:
                z_data[i] = str(z)
                print(z_data[i])
    
        np.savetxt('redshift.txt', np.transpose([id, z_data]) , fmt='%s', delimiter='\t')
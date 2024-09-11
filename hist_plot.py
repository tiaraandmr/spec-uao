#!/usr/bin/python python3

import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt
plt.style.use('presentation.mplstyle')
from astropy.io import fits

# Define the home directory
home = '/Users/saraswati/Documents/Work/spec-uao/'

# Read the object list
slits_path = home+'slits_reduced_z.fits'
slits = Table.read(slits_path)

# Define the arrays
z = slits['SPEC_Z']
redshift = z[~np.isnan(z)]
has_redshift = ~np.isnan(z)
no_redshift = np.isnan(z)
mag_i = slits['MAG_I']

# Plot the single histogram
plt.hist(mag_i, bins=6, range=(21,24), color='black',histtype='step',lw=2)
plt.hist(mag_i[has_redshift], bins=6, range=(21,24), color='crimson',histtype='step',lw=2)
plt.hist(mag_i[no_redshift], bins=6, range=(21,24), color='teal',histtype='step',lw=2)
plt.ylabel(r'$\#$ Spectra')
plt.xlabel(r'HSC $i$ (AB)')

# Add the legend
legend_1 = plt.plot([],[], color='crimson', linestyle='-', label=r'$z$ recovered (56)')
legend_2 = plt.plot([],[], color='teal', linestyle='-', label=r'No-$z$ (80)')

# Plot the legend only
plt.legend(prop={'size': 13})
plt.savefig('hist.pdf', bbox_inches="tight")
plt.close()

# Subplots parameters
fig,axes = plt.subplots(2,2,figsize=(12,8),sharex='col',sharey='row',gridspec_kw={'height_ratios': [1,3],'width_ratios':[3,1]})
fig.subplots_adjust(hspace=0,wspace=0)

# Plot the bottom left scatter plot
axes[1,0].scatter(z, mag_i, color='black', s=100)
axes[1,0].set_xlabel(r'MMT $z_{\mathrm{spec}}$')
axes[1,0].set_ylabel(r'HSC $i$ (AB)')
axes[1,0].set_ylim(20.5,24)

# Plot the top left histogram
axes[0,0].hist(redshift,bins=8, range=(0,2), color='black',histtype='step',lw=2)
axes[0,0].set_ylabel(r'$\#$ Spectra', labelpad=20)
axes[0,0].set_xlabel(r'MMT $z_{\mathrm{spec}}$', labelpad=20)
axes[0,0].xaxis.set_tick_params(labelbottom=True)
axes[0,0].xaxis.set_label_position("top")
axes[0,0].xaxis.tick_top()

# Plot the bottom right histogram
axes[1,1].hist(mag_i,bins=6, range=(21,24), color='black',histtype='step',lw=2, orientation='horizontal')
axes[1,1].set_ylabel(r'HSC $i$ (AB)', rotation=270, labelpad=25)
axes[1,1].set_xlabel(r'$\#$ Spectra')
axes[1,1].yaxis.set_label_position("right")
axes[1,1].yaxis.set_tick_params(labelbottom=True)
axes[1,1].yaxis.tick_right()

# Turn off the top right plot
axes[0,1].axis('off')

# Save the figure
plt.savefig('hist_all.pdf', bbox_inches="tight")
plt.close()




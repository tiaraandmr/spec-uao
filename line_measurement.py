#!/usr/bin/env python3

# Import packages
import numpy as np
from astropy.table import Table
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

import matplotlib.pyplot as plt
from astropy.io import fits
import argparse

# Define parser
parser = argparse.ArgumentParser()
parser.add_argument('--id', default='41231836365464984')
args = parser.parse_args()
results_id = args.id

# Define path to .FITS file
results_path = f'/Users/saraswati/Documents/Work/spec-uao/results/{results_id}-results.fits'
results = fits.open(results_path)

#.FITS data
data = results[2].data

#.FITS header
header = results[2].header

# Line flux value
ne_III = np.average(data['AGN_[NeIII]_3869.86_Flux'])
o_II = np.average(data['SF_[OII]_3728.48_Flux'])
o_III = np.average(data['AGN_[OIII]_5008.24_Flux'])
h_beta = np.average(data['Balmer_HI_4862.68_Flux'])  
#n_II = np.average(data['AGN_[NII]_6585.27_Flux'])
#h_alpha = np.average(data['Balmer_HI_6564.61_Flux'])

# Line width value
o_III_w = np.average(data['AGN_[OIII]_5008.24_REW'])
#h_alpha_w = np.average(data['Balmer_HI_6564.61_REW'])

# Line dispersion value
o_III_d = np.average(data['AGN_[OIII]_5008.24_Dispersion'])
#h_alpha_d = np.average(data['Balmer_HI_6564.61_Dispersion'])

ratio_1 = np.log10(ne_III/o_II)
ratio_2 = np.log10(o_III/h_beta)
#ratio_3 = np.log10(n_II/h_alpha)

print(results_id)
print(r'LOG([NeIII]/[OII]) = ', ratio_1)
print(r'LOG([OIII]/H-beta) = ', ratio_2)
#print(r'LOG([NII]/H-alpha) = ', ratio_3)
print(r'[OIII] REW = ', o_III_w)
#print(r'H_alpha REW = ', h_alpha_w)
print(r'[OIII] Dispersion = ', o_III_d)
#print(r'H_alpha Dispersion = ', h_alpha_d)

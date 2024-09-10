#!/usr/bin/python python3

# Import packages
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
from astropy.table import Table
from spectres import spectres
from astropy.convolution import convolve

# Set the font style to LaTeX
rc('font', **{'family': 'serif', 'serif': ['lmodern']})
rc('text', usetex=True)

# Define home directory
home = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/'

# Load the flux sampling template
sampling_path = '/Users/saraswati/Documents/Work/spec-uao/calibration_star/BD+33d2642_20230708_sampling.fits'
sampling = Table.read(sampling_path)

# Load the observed spectrum
spec_obs_path = '/Users/saraswati/Documents/Work/spec-uao/clean_spec/41231853545355348.fits'
spec_obs = Table.read(spec_obs_path)

# create a new wavelength array for both spectra
new_wav = spec_obs['WAVELENGTH']
flux = spec_obs['FLUX']

# resample the 'sampling' spectra based on the new wavelength array
flux_sampling_resample = spectres(new_wav, sampling['WAVELENGTH'], sampling['FLUX'], verbose=False)

# Apply the flux sampling to the observed spectra
flux_calibrated = spec_obs['FLUX'] * flux_sampling_resample

# automatically set the x-limit for every spectra
good = np.invert(np.isnan(flux_calibrated))

#define line species
species = np.array([['[SII]', 6716.435709204644], ['[SII]', 6730.811841618589], ['[NII]', 6583.451474377235], ['[NII]', 6548.050991532048], 
                    ['[OIII]', 5006.843330296288], ['[OIII]', 4958.91106815638], ['[OIII]', 4363.2096874449635], ['[NeIII]', 3868.7632032493057], 
                    ['[NeV]', 3425.8676240424247], ['[NeV]', 3345.828076914398], ['[OI]', 6300.303551730963], ['[OI]', 6363.776506749638], 
                    ['[OII]', 3727.4199178007175], ['HI', 6562.797027356974], ['HI', 4861.321979760415], ['HI', 4340.459677187083], 
                    ['CIV', 1548.8576020874154], ['HeII', 1639.7908206670402], ['CIII]', 1908.1020362243723], ['MgII', 2798.2921038031536]])

slits = Table.read('/Users/saraswati/Documents/Work/spec-uao/slits_reduced.fits')
object_ID = slits['ID'].astype('str')
#z = slits['BEST_Z']
z = 0.285

#smooth the spectra
N = 5
spec_convolve = convolve(flux_calibrated, np.ones(N), boundary='extend', normalize_kernel=True)
spec_convolve_orig = convolve(spec_obs['FLUX'], np.ones(N), boundary='extend', normalize_kernel=True)

#define the figure and font size
plt.rcParams.update({'font.size': 17})
fig, ax = plt.subplots(figsize=(17,5))

#plot the spectrum and set the x-limit
ax.semilogy(new_wav, spec_convolve, color='firebrick', linewidth=2.0, drawstyle='steps-mid', label='Calibrated \& Smoothed (N=5)', alpha=0.7)
ax.semilogy(new_wav, spec_convolve_orig, color='steelblue', linewidth=2.0, drawstyle='steps-mid', label='Not Calibrated, Smoothed (N=5)', alpha=0.7)
ax.legend(bbox_to_anchor=[0.48,0.3])
    
#define the approximate redshift (z)
best_z = z

#check if lines are in the range of the observed spectrum
for j in range(len(species)):
    data = species[j][1].astype('float64')
    x = data*(1+best_z)

    if ((x < new_wav[-1]) and (x > new_wav[0])):
        ax.vlines(x, ymin=flux.min(), ymax=flux.max(), color = 'black', ls = '--') 
        ax.text(x, 0.3*10*flux.max(), species[j][0], rotation=90, fontsize=14, ha='center', va='center')

#setting the x- and y-axis limit
ax.set_xlim(new_wav[good].min(), new_wav[good].max())
ax.set_ylim(flux.min(), 100*flux.max())

#plot a rest frame wavelength based on the best_z
ax2 = ax.secondary_xaxis('top',functions=(lambda lam: lam/(1+best_z), lambda lam: lam*(1+best_z))) 
ax2.set_xlabel(r'Rest Wavelength [$ \rm \AA$]', labelpad=10)
    
#define x and y label and plot title
ax.set_xlabel(r'Observed Wavelength [$ \rm \AA$]', labelpad=10)
ax.set_ylabel(r'Flux [$\mathrm{10^{-17}\ erg\ cm^{-2}\ s^{-1}\ \AA^{-1}}$]', labelpad=10)
title = '41231853545355348'
ax.set_title(title+r', $z = {:.3f}$'.format(best_z), pad=20)

#save the plot as .pdf
plt.savefig(home+title+'_calibrated_smooth_n5.pdf', bbox_inches="tight")
plt.close()
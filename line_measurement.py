#!/usr/bin/env python3

# Import packages
import os
import numpy as np
import pandas as pd
from astropy.table import Table
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

import matplotlib.pyplot as plt
from astropy.io import fits

# Open results file
home = '/Users/saraswati/Documents/Work/spec-uao/'
results_path = home+'results/GELATO-results.fits'
results = Table.read(results_path, 1)

ratio_1 = pd.Series([], dtype='float64')
ratio_1_err = pd.Series([], dtype='float64')
ratio_2 = pd.Series([], dtype='float64')
ratio_2_err = pd.Series([], dtype='float64')
ratio_3 = pd.Series([], dtype='float64')
ratio_3_err = pd.Series([], dtype='float64')
name_clean = pd.Series([], dtype='int')
rest_color_val = pd.Series([], dtype='float64')
rest_color_err = pd.Series([], dtype='float64')

# Open rest frame color file
rest_color_path = home+'SED-concat.fits'
rest_color = Table.read(rest_color_path,1)

for i in range(len(results)):

    # Line flux value
    ne_III = results['AGN_[NeIII]_3869.86_Flux'][i]
    ne_III_err = results['AGN_[NeIII]_3869.86_Flux_err'][i]
    o_II = results['SF_[OII]_3728.48_Flux'][i]
    o_II_err = results['SF_[OII]_3728.48_Flux_err'][i]
    o_III = results['AGN_[OIII]_5008.24_Flux'][i]
    o_III_err = results['AGN_[OIII]_5008.24_Flux_err'][i]
    h_beta = results['Balmer_HI_4862.68_Flux'][i] 
    h_beta_err = results['Balmer_HI_4862.68_Flux_err'][i]
    n_II = results['AGN_[NII]_6585.27_Flux'][i]
    n_II_err = results['AGN_[NII]_6585.27_Flux_err'][i]
    h_alpha = results['Balmer_HI_6564.61_Flux'][i]
    h_alpha_err = results['Balmer_HI_6564.61_Flux_err'][i]

    # List the line
    lines_1 = [ne_III, o_II] 
    lines_2 = [o_III, h_beta]
    lines_3 = [n_II, h_alpha]

    for j in range(len(lines_1)):
        if lines_1[0] != 'nan' and lines_1[1] != 'nan':
            if lines_1[0] > 0 and lines_1[1] > 0:
                ratio_1[i] = np.log10(ne_III/o_II)
                ratio_1_err[i] = 0.434 * np.sqrt((ne_III_err/ne_III)+(o_II_err/o_II))
            else:
                ratio_1[i] = np.nan
                ratio_1_err[i] = np.nan
        else:
            ratio_1[i] = np.nan
            ratio_1_err[i] = np.nan

    for k in range(len(lines_2)):
        if lines_2[0] != 'nan' and lines_2[1] != 'nan':
            if lines_2[0] > 0 and lines_2[1] > 0:
                ratio_2[i] = np.log10(o_III/h_beta)
                ratio_2_err[i] = 0.434 * np.sqrt((o_III_err/o_III)+(h_beta_err/h_beta))
            else:
                ratio_2[i] = np.nan
                ratio_2_err[i] = np.nan
        else:
            ratio_2[i] = np.nan
            ratio_2_err[i] = np.nan
    
    for l in range(len(lines_3)):
        if lines_3[0] != 'nan' and lines_3[1] != 'nan':
            if lines_3[0] > 0 and lines_3[1] > 0:
                ratio_3[i] = np.log10(n_II/h_alpha)
                ratio_3_err[i] = 0.434 * np.sqrt((n_II_err/n_II)+(h_alpha_err/h_alpha))
            else:
                ratio_3[i] = np.nan
                ratio_3_err[i] = np.nan
        else:
            ratio_3[i] = np.nan
            ratio_3_err[i] = np.nan
    
    name_clean[i] = results['Name'][i].rsplit('.', 1)[0]

    for m in range(len(rest_color)):
        if int(rest_color['ID'][m]) == int(name_clean[i]):
            rc_z = rest_color['rest_sdss_z'][m]
            err_rc_z = rest_color['rest_sdss_z_err'][m]
            rc_g = rest_color['rest_sdss_g'][m]
            err_rc_g = err_rc_z = rest_color['rest_sdss_g_err'][m]
            rest_color_val[i] = 2.5 * np.log10(rc_z/rc_g)
            rest_color_err[i] = 2.5 * 0.434 * np.sqrt((err_rc_z/rc_z)+(err_rc_g/rc_g))
        else:
            pass

#save value into a DataFrame
ratio = pd.DataFrame()
ratio.insert(0, "ID", name_clean)
ratio.insert(1, "Ne_III/O_II_Flux_Ratio", ratio_1)
ratio.insert(2, "Ne_III/O_II_Flux_Ratio_Err", ratio_1_err)
ratio.insert(3, "O_III/H_Beta_Flux_Ratio", ratio_2)
ratio.insert(4, "O_III/H_Beta_Flux_Ratio_Err", ratio_2_err)
ratio.insert(5, "N_II/H_Alpha_Flux_Ratio", ratio_3)
ratio.insert(6, "N_II/H_Alpha_Flux_Ratio_Err", ratio_3_err)
ratio.insert(7, "Rest_SDSS_Color", rest_color_val)
ratio.insert(8, "Rest_SDSS_Color_Err", rest_color_err)

#ratio_tab = Table([name_clean, ratio_1, ratio_1_err, ratio_2, ratio_2_err,
#                   ratio_3, ratio_3_err, rest_color_val, rest_color_err], 
#                   names=('ID', 'Ne_III/O_II_Flux_Ratio', 'Ne_III/O_II_Flux_Ratio_Err', 
#                          'O_III/H_Beta_Flux_Ratio', 'O_III/H_Beta_Flux_Ratio_Err', 'N_II/H_Alpha_Flux_Ratio',
#                          'N_II/H_Alpha_Flux_Ratio_Err', 'Rest_SDSS_Color', 'Rest_SDSS_Color_Err'))

#convert DataFrame into Table
ratio_tab = Table.from_pandas(ratio)

os.chdir(home+'results')
ratio_tab.write('line_ratio.fits', overwrite=True)

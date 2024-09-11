#!/usr/bin/python python3

import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.style.use('presentation.mplstyle')
from astropy.io import fits

# Define home directory
home = '/home/tiara/Documents/Work/spec-uao/'

# Define results directory
results = home+'analysis_plot/'

# Load the line ratio table
line_ratio = home+'line_ratio_broad.fits'
ratio = Table.read(line_ratio, 1)

# Define the empirical lines
emp_x = np.log10(np.linspace(0.01,4,1000))
emp_line_tbt = -1.2 * emp_x - 0.4
emp_line_bpt = 0.61 / (emp_x - 0.05) + 1.3

# Plot the TBT diagram (single plot)
plt.fill_between(emp_x, -0.5, emp_line_tbt, facecolor='deepskyblue', alpha=0.2)
plt.fill_between(emp_x, emp_line_tbt, 2, facecolor='crimson', alpha=0.2)

for i in range(len(ratio)):
    if ratio['P_AGN'][i] > 80:
        #axes[0,0].plot(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_AGN'][i], linestyle='None', color='crimson', marker='.', markersize=10, label='AGN Color')
        plt.errorbar(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_AGN'][i], xerr=float(ratio['Ne_III/O_II_Flux_Ratio_Err'][i]), yerr=float(ratio['Rest_SDSS_Color_Err_AGN'][i]),
                          color='black', marker='.', markersize=10, capsize=2, ecolor='black', elinewidth=1)
    else:
        #axes[0,0].plot(ratio['Ne_III/O_II_Flux_Ratio'], ratio['Rest_SDSS_Color_Gal'], linestyle='None', color='deepskyblue', marker='.', markersize=10, label='SF Color')
        plt.errorbar(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_Gal'][i], xerr=float(ratio['Ne_III/O_II_Flux_Ratio_Err'][i]), yerr=float(ratio['Rest_SDSS_Color_Err_Gal'][i]),
                          color='black', marker='.', markersize=10, capsize=2, ecolor='black', elinewidth=1)

plt.plot(emp_x, emp_line_tbt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
plt.ylabel(r'$^{0.0}(g-z)$')
plt.xlabel(r'log$_{10}$([NeIII]/[OII])')
plt.xlim(-2,0.3)
plt.ylim(-0.5,2)
plt.title('TBT ($z<1.2$)')

plt.savefig(results+'tbt_1.pdf', bbox_inches="tight")
plt.close()

# Plot the PBT diagram (single plot)
for i in range(len(ratio)):
    if ratio['N_II/H_Alpha_Flux_Ratio'][i] != np.nan and ratio['O_III/H_Beta_Flux_Ratio'][i] != np.nan:
        test_value = 0.61 / (float(ratio['N_II/H_Alpha_Flux_Ratio'][i]) - 0.05) + 1.3
        if test_value > ratio['O_III/H_Beta_Flux_Ratio'][i]:
            plt.plot(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], color='black', linestyle='None', marker='.', markersize=10)
            plt.errorbar(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], xerr=ratio['N_II/H_Alpha_Flux_Ratio_Err'][i], yerr=ratio['O_III/H_Beta_Flux_Ratio_Err'][i],
                          color='black', capsize=2, ecolor='black', elinewidth=1)
        else:
            plt.plot(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], color='black', linestyle='None', marker='.', markersize=10)
            plt.errorbar(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], xerr=ratio['N_II/H_Alpha_Flux_Ratio_Err'][i], yerr=ratio['O_III/H_Beta_Flux_Ratio_Err'][i],
                          color='black', capsize=2, ecolor='black', elinewidth=1)
    else:
        pass

plt.plot(emp_x, emp_line_bpt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
plt.fill_between(emp_x, -1, emp_line_bpt, facecolor='deepskyblue', alpha=0.2)
plt.fill_between(emp_x, emp_line_bpt, 1.5, facecolor='crimson', alpha=0.2)
plt.ylabel(r'log$_{10}$([OIII]/H$\beta$)')
plt.xlabel(r'log$_{10}$([NII]/H$\alpha$)')
plt.xlim(-2,0)
plt.ylim(-1,1.5)
plt.title('BPT ($z<1.2$)')

plt.savefig(results+'pbt_1.pdf', bbox_inches="tight")
plt.close()

# Plot the PBT, TBT, and Ne V with SNR > 3 (4-panel plots)
plt.rcParams.update({'font.size': 17})
fig, axes = plt.subplots(2, 2, figsize=(12,10))
fig.subplots_adjust(wspace=0.25, hspace=0.4)

color = []

for i in range(len(ratio)):
    if ratio['P_AGN'][i] > 80:
        color_code = 'crimson'
        color.append(color_code)
    else:
        color_code = 'deepskyblue'
        color.append(color_code)

axes[0,1].fill_between(emp_x, -0.5, emp_line_tbt, facecolor='deepskyblue', alpha=0.2)
axes[0,1].fill_between(emp_x, emp_line_tbt, 2, facecolor='crimson', alpha=0.2)

for a, b, c in zip(ratio['Ne_III/O_II_Flux_Ratio'], ratio['Rest_SDSS_Color_AGN'], 
                         ratio['Rest_SDSS_Color_Gal']):
    axes[0,1].scatter(float(a), float(b), color='crimson', linestyle='None', marker='x', linewidths=2, s=80)
    axes[0,1].scatter(float(a), float(c), color='deepskyblue', linestyle='None', marker='x', linewidths=2, s=80)
    axes[0,1].plot([a,a], [b,c], color='black', linestyle='-', alpha=0.7)

axes[0,1].plot(emp_x, emp_line_tbt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
axes[0,1].set_ylabel(r'$^{0.0}(g-z)$')
axes[0,1].set_xlabel(r'log$_{10}$([NeIII]/[OII])')
axes[0,1].set_xlim(-2,0.3)
axes[0,1].set_ylim(-0.5,2)
axes[0,1].set_title('TBT ($z<1.2$)')

legend_1 = axes[0,1].scatter([],[], color='crimson', linestyle='None', marker='x', linewidths=2, s=80, label='AGN Color')
legend_2 = axes[0,1].scatter([],[], color='deepskyblue', linestyle='None', marker='x', linewidths=2, s=80, label='SF Color')
axes[0,1].legend(loc=4, prop={'size': 13}, ncol=1)

for i in range(len(ratio)):
    if ratio['P_AGN'][i] > 80:
        #axes[0,0].plot(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_AGN'][i], linestyle='None', color='crimson', marker='.', markersize=10, label='AGN Color')
        axes[0,0].errorbar(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_AGN'][i], xerr=float(ratio['Ne_III/O_II_Flux_Ratio_Err'][i]), yerr=float(ratio['Rest_SDSS_Color_Err_AGN'][i]),
                          color='black', marker='.', markersize=10, capsize=2, ecolor='black', elinewidth=1)
    else:
        #axes[0,0].plot(ratio['Ne_III/O_II_Flux_Ratio'], ratio['Rest_SDSS_Color_Gal'], linestyle='None', color='deepskyblue', marker='.', markersize=10, label='SF Color')
        axes[0,0].errorbar(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_Gal'][i], xerr=float(ratio['Ne_III/O_II_Flux_Ratio_Err'][i]), yerr=float(ratio['Rest_SDSS_Color_Err_Gal'][i]),
                          color='black', marker='.', markersize=10, capsize=2, ecolor='black', elinewidth=1)

agn_tbt = []

for i in range(len(ratio)):
    if ratio['Ne_III/O_II_Flux_Ratio'][i] != np.nan:
        test_value = -1.2 * ratio['Ne_III/O_II_Flux_Ratio'][i] - 0.4
        if ratio['P_AGN'][i] > 80:
            if test_value < ratio['Rest_SDSS_Color_AGN'][i]:
                id_val = ratio['ID'][i]
                agn_tbt.append(id_val)
        else:
            if test_value < ratio['Rest_SDSS_Color_Gal'][i]:
                id_val = ratio['ID'][i]
                agn_tbt.append(id_val)
    
for i in range(len(agn_tbt)):
    for j in range(len(ratio)):
        if agn_tbt[i] == ratio['ID'][j]:
            if ratio['Mg_II_Dispersion'][j] > 500:
                print(ratio['ID'][j])

axes[0,0].fill_between(emp_x, -0.5, emp_line_tbt, facecolor='deepskyblue', alpha=0.2)
axes[0,0].fill_between(emp_x, emp_line_tbt, 2, facecolor='crimson', alpha=0.2)

axes[0,0].plot(emp_x, emp_line_tbt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
axes[0,0].set_ylabel(r'$^{0.0}(g-z)$')
axes[0,0].set_xlabel(r'log$_{10}$([NeIII]/[OII])')
axes[0,0].set_xlim(-2,0.3)
axes[0,0].set_ylim(-0.5,2)
axes[0,0].set_title('TBT ($z<1.2$)')

#legend_1 = axes[0,0].plot([],[], color='crimson', linestyle='None', marker='.', markersize=10, label=r'P(AGN) $>$ 80')
#legend_2 = axes[0,0].plot([],[], color='deepskyblue', linestyle='None', marker='.', markersize=10, label=r'P(AGN) $<$ 80')
#axes[0,0].legend(loc=4, prop={'size': 13}, ncol=1)
#axes[0,0].legend(loc=8, prop={'size': 13}, ncol=2)

for i in range(len(ratio)):
    if ratio['SNR_Ne_V'][i] > 3:
        axes[1,0].plot(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_AGN'][i], color='black', linestyle='None', marker='.', markersize=10)
        axes[1,0].errorbar(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_AGN'][i], xerr=float(ratio['Ne_III/O_II_Flux_Ratio_Err'][i]), yerr=float(ratio['Rest_SDSS_Color_Err_AGN'][i]),
                          color='black', capsize=2, ecolor='black', elinewidth=1)

axes[1,0].fill_between(emp_x, -0.5, emp_line_tbt, facecolor='deepskyblue', alpha=0.2)
axes[1,0].fill_between(emp_x, emp_line_tbt, 2, facecolor='crimson', alpha=0.2)

axes[1,0].plot(emp_x, emp_line_tbt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
axes[1,0].set_ylabel(r'$^{0.0}(g-z)$')
axes[1,0].set_xlabel(r'log$_{10}$([NeIII]/[OII])')
axes[1,0].set_xlim(-2,0.3)
axes[1,0].set_ylim(-0.5,2)
axes[1,0].set_title('[NeV] SNR $> 3$')

for i in range(len(ratio)):
    if ratio['N_II/H_Alpha_Flux_Ratio'][i] != np.nan and ratio['O_III/H_Beta_Flux_Ratio'][i] != np.nan:
        test_value = 0.61 / (float(ratio['N_II/H_Alpha_Flux_Ratio'][i]) - 0.05) + 1.3
        if test_value > ratio['O_III/H_Beta_Flux_Ratio'][i]:
            axes[1,1].plot(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], color='black', linestyle='None', marker='.', markersize=10)
            axes[1,1].errorbar(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], xerr=ratio['N_II/H_Alpha_Flux_Ratio_Err'][i], yerr=ratio['O_III/H_Beta_Flux_Ratio_Err'][i],
                          color='black', capsize=2, ecolor='black', elinewidth=1)
        else:
            axes[1,1].plot(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], color='black', linestyle='None', marker='.', markersize=10)
            axes[1,1].errorbar(ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['O_III/H_Beta_Flux_Ratio'][i], xerr=ratio['N_II/H_Alpha_Flux_Ratio_Err'][i], yerr=ratio['O_III/H_Beta_Flux_Ratio_Err'][i],
                          color='black', capsize=2, ecolor='black', elinewidth=1)
    else:
        pass

axes[1,1].plot(emp_x, emp_line_bpt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
axes[1,1].fill_between(emp_x, -1, emp_line_bpt, facecolor='deepskyblue', alpha=0.2)
axes[1,1].fill_between(emp_x, emp_line_bpt, 1.5, facecolor='crimson', alpha=0.2)
axes[1,1].set_ylabel(r'log$_{10}$([OIII]/H$\beta$)')
axes[1,1].set_xlabel(r'log$_{10}$([NII]/H$\alpha$)')
axes[1,1].set_xlim(-2,0)
axes[1,1].set_ylim(-1,1.5)
axes[1,1].set_title('BPT ($z<1.2$)')


fig.savefig(results+'tbt_all_color_ver.pdf', bbox_inches="tight")
fig.clf()
plt.close(fig)
#!/usr/bin/python python3

import numpy as np
import numpy.ma as ma
import pandas as pd
from astropy.table import Table

import matplotlib.pyplot as plt
plt.style.use('presentation.mplstyle')
from astropy.io import fits

home = '/Users/saraswati/Documents/Work/spec-uao/'
line_ratio = home+'line_ratio_broad.fits'
ratio = Table.read(line_ratio, 1)

emp_x = np.log10(np.linspace(0.01,4,1000))
emp_line_tbt = -1.2 * emp_x - 0.4
emp_line_bpt = 0.61 / (emp_x - 0.05) + 1.3

#define the figure and font size
#plt.figure(figsize=(7,5))
plt.rcParams.update({'font.size': 17})
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,5))

#color = []
#    if ratio['P_AGN'][i] > 80:
#        color_code = 'crimson'
#        color.append(color_code)
#    else:
#        color_code = 'deepskyblue'
#        color.append(color_code)

#for a, b, c, d, e in zip(ratio['Ne_III/O_II_Flux_Ratio'], ratio['Rest_SDSS_Color_AGN'], 
#                         ratio['Rest_SDSS_Color_Gal'], marker_agn, marker_gal):
#    plt.plot([a,a], [b,c], color='black', linestyle='-')
#for a, b, c, d, e in zip(ratio['Ne_III/O_II_Flux_Ratio'], ratio['Rest_SDSS_Color'], ratio['Ne_III/O_II_Flux_Ratio_Err'], 
#                         ratio['Rest_SDSS_Color_Err'], color):
#    plt.plot(a, b, color=e, linestyle='None', marker='.', markersize=12)
#    plt.errorbar(a, b, xerr=float(c), yerr=float(d), color=e, capsize=2, ecolor='black', elinewidth=0.3)

ax1.plot([ratio['Ne_III/O_II_Flux_Ratio'], ratio['Ne_III/O_II_Flux_Ratio']], [ratio['Rest_SDSS_Color_AGN'], ratio['Rest_SDSS_Color_Gal']],
         'k-', alpha=0.7)
ax1.plot(ratio['Ne_III/O_II_Flux_Ratio'], ratio['Rest_SDSS_Color_AGN'], color='crimson', linestyle='None', marker='x', markersize=8, mew=3, label='AGN Color')
ax1.plot(ratio['Ne_III/O_II_Flux_Ratio'], ratio['Rest_SDSS_Color_Gal'], color='deepskyblue', linestyle='None', marker='.', markersize=12, label='SF Color')

ax1.fill_between(emp_x, -0.5, emp_line_tbt, facecolor='deepskyblue', alpha=0.2)
ax1.fill_between(emp_x, emp_line_tbt, 2, facecolor='crimson', alpha=0.2)

ax1.plot(emp_x, emp_line_tbt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
ax1.set_ylabel(r'$^{0.0}(g-z)$')
ax1.set_xlabel(r'log$_{10}$([NeIII]/[OII])')
ax1.set_xlim(-2,0.3)
ax1.set_ylim(-0.5,2)
ax1.set_title('TBT ($z<1.2$)')
#ax1.legend(loc=8, prop={'size': 13}, ncol=2)

for i in range(len(ratio)):
    if ratio['SNR_Ne_V'][i] > 3:
        ax2.plot([ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Ne_III/O_II_Flux_Ratio'][i]], [ratio['Rest_SDSS_Color_AGN'][i], ratio['Rest_SDSS_Color_Gal'][i]],
         'k-', alpha=0.7)
        ax2.plot(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_AGN'][i], color='crimson', linestyle='None', marker='x', markersize=8, mew=3, label='AGN Color')
        ax2.plot(ratio['Ne_III/O_II_Flux_Ratio'][i], ratio['Rest_SDSS_Color_Gal'][i], color='deepskyblue', linestyle='None', marker='.', markersize=12, label='SF Color')

ax2.fill_between(emp_x, -0.5, emp_line_tbt, facecolor='deepskyblue', alpha=0.2)
ax2.fill_between(emp_x, emp_line_tbt, 2, facecolor='crimson', alpha=0.2)

ax2.plot(emp_x, emp_line_tbt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
ax2.set_ylabel(r'$^{0.0}(g-z)$')
ax2.set_xlabel(r'log$_{10}$([NeIII]/[OII])')
ax2.set_xlim(-2,0.3)
ax2.set_ylim(-0.5,2)
ax2.set_title('[NeV] SNR $> 3$')
#ax2.legend(loc=8, prop={'size': 13}, ncol=2)

#ax2.plot([ratio['N_II/H_Alpha_Flux_Ratio'][i], ratio['N_II/H_Alpha_Flux_Ratio'][i]], [ratio['Rest_SDSS_Color_AGN'][i], ratio['Rest_SDSS_Color_Gal'][i]],
#         'k-', alpha=0.7)
ax3.plot(ratio['N_II/H_Alpha_Flux_Ratio'], ratio['O_III/H_Beta_Flux_Ratio'], color='crimson', linestyle='None', marker='x', markersize=8, mew=3, label='AGN Color')
ax3.plot(ratio['N_II/H_Alpha_Flux_Ratio'], ratio['O_III/H_Beta_Flux_Ratio'], color='deepskyblue', linestyle='None', marker='.', markersize=12, label='SF Color')

ax3.plot(emp_x, emp_line_bpt, color='black', linewidth=2.0, linestyle='--', alpha=0.5)
ax3.fill_between(emp_x, -1, emp_line_bpt, facecolor='deepskyblue', alpha=0.2)
ax3.fill_between(emp_x, emp_line_bpt, 1.5, facecolor='crimson', alpha=0.2)
ax3.set_ylabel(r'log$_{10}$([OIII]/H$\beta$)')
ax3.set_xlabel(r'log$_{10}$([NII]/H$\alpha$)')
ax3.set_xlim(-2,0)
ax3.set_ylim(-1,1.5)
ax3.set_title('BPT ($z<1.2$)')


fig.savefig('tbt_all.pdf', bbox_inches="tight")
plt.close(fig)
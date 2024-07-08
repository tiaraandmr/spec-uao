import os, glob, shutil, subprocess as sub
from pathlib import Path

import numpy as np
from astropy.table import Table
from astropy.convolution import convolve

from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)
import matplotlib.pyplot as plt
#from matplotlib.widgets import TextBox
from astropy.io import fits

#def visualizeGraph(expr):
#    ydata = eval(expr)
#    spec_convolve.set_ydata(ydata)
#    ax.relim()
#    ax.autoscale_view()
#    plt.draw()
 
# Adding TextBox to graph
#graphBox = fig.add_axes([0.1, 0.05, 0.8, 0.075])
#txtBox = TextBox(graphBox, r'Plot: ')
#txtBox.on_submit(visualizeGraph)
#txtBox.set_val('0.81')
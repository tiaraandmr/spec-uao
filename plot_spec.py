#!/usr/bin/python python3

import os, subprocess as sub
from pathlib import Path
import numpy as np
import numpy.ma as ma
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['lmodern']})
rc('text', usetex=True)

from spectres import spectres
from specutils import Spectrum1D
from astropy.io import fits
from astropy.wcs import WCS
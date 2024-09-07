# Data Calibration and Analysis Code for MMT Observatory Binospec Spectra

## Dependencies
To properly run each scripts in this repository, the following packages are required:
- NumPy
- Astropy
- Pandas
- Spectres

To properly generate the plots, matplotlib must be installed, and a version of LaTeX must be available as well. LaTeX can be installed following the instructions [here](https://www.latex-project.org/get/).

It is advised to clone the entire spec-uao repository. This is done by using the following commands to clone the directory over HTTPS.

```
cd /path/to/directory
git clone https://github.com/tiaraandmr/spec-uao.git
```

Then the dependencies can be installed. It is recommended to install conda (or another wrapper). A conda environment with nearly all the dependencies can be installed via the provided "environment.yml" file. If you do not wish to use conda, you can install the dependencies enumerated in the "environment.yml" file.

```
cd /path/to/spec-uao/directory
conda env create -f environment.yml
```

The environment will then be installed under the name spec-uao and can then be activated.

```
conda activate spec-uao
```

## Running the Scripts
### Spectral Cleaning
The spectra is reduced with the Binospec pipeline described [here](https://bitbucket.org/chil_sai/binospec/wiki/Home). The 1D spectra is now split into individual fits file, one spectrum per file. They are distributed bundled up into tar files, whose contents are the 1D files with file names derived from object names submitted with the mask design. For the example here, obj_abs_1D.tar folder is the default folder used. The spectra from the pipeline still need further cleaning from the sky lines, bad data, and calibration artifact. This script produce new fits table with the cleaned spectrum for each file (clean_spec) and preliminary plots (plot_spec) for a quicklook. After the object ID is obtained, this script also cross-match them with the bigger object catalog "slits.fits", and then produced new catalog with only the object in the sample "slits_reduced.fits".

```
python3 clean_spec.py
```

It is advised to always check the path to the fits file written on the beginning of every scripts.

### Spectral Smoothing
To make it easier to identify the lines, the clean spectra are smoothed to a specific degree. A preliminary redshift (z) values will be plotted together with the rest wavelength of the emission lines (smooth_spec). 

```
python3 smooth_spec.py
```
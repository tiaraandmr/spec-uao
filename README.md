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

It is advised to always check the path to the fits file written on the beginning of every scripts. If there are observation from two different nights for one objects, this script will also combined them with a weighted average. 

### Spectral Smoothing
To make it easier to identify the lines, the clean spectra are smoothed to a specific degree using [Spectres](https://github.com/ACCarnall/spectres?tab=readme-ov-file). A preliminary redshift (z) values will be plotted together with the rest wavelength of the emission lines (smooth_spec). 

```
python3 smooth_spec.py
```

### Flux Calibration
After the spectra is cleaned, the flux needs to be calibrated using a F-star or known spectrophotometric standard star. In this example, we use spectrophotometric standard star with observation date closest to each field. The spectra of the standard star needs to be resampled to match the wavelength range of the object spectra. This scripts will resample the standard star (input file is a .txt) then produce a fits file ready for calibration and a quicklook plot. All the results will be generated inside the calibration_star folder.

```
python3 flux_sampling.py
```
The next step is applying the calibration to all of the spectra. This scripts will apply the flux calibration and then save the results in the calibrated_spec folder as fits file. Aside from that, fits file for fitting purpose with [GELATO](https://github.com/TheSkyentist/GELATO) will also be generated.

```
python3 flux_calib_all.py
```

### Redshift Estimation
We can determine the redshift by searching for a prominent emission or absorption lines. To improve the visual inspection, we make use of [redshifting](https://github.com/sdjohnson-astro/redshifting) code which performs a grid-search using eigenspectrum templates from [Bolton et al. (2012)](https://iopscience.iop.org/article/10.1088/0004-6256/144/5/144) to determine redshift from optical spectroscopy. The example of parameters obtained with redshifting are given in plot_redshifting and redshifting_redshift.fits. The following script open an interactive plots of the input spectra and then record the last redshift input by the user to a file, along with the object ID (redshift.txt). Object ID has to be specified when running the script, for example:

```
python3 interactive_spec.py --id 42300690516701462
```
After determining the redshift, [GELATO](https://github.com/TheSkyentist/GELATO) will be used to fit the spectra. An example of GELATO results are given in results_high_z_broad (for spectra with z > 1.2 and broad lines parameter), results_low_z_broad (for spectra with z < 1.2 and broad lines parameter), and results_low_z_narrow (for spectra with z < 1.2 and narrow lines parameter).

### Line Measurements
This script will get the parameters needed to create an AGN diagnostic diagram, which is [PBT](https://iopscience.iop.org/article/10.1086/130766) and [TBT](https://iopscience.iop.org/article/10.1088/0004-637X/742/1/46) diagram, from the fitting results. After line flux and dispersion were obtained, line ratio will be calculated and stored into a fits file (line_ratio_broad.fits).

```
python3 line_measurement.py
```


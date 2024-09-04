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
### Cleaning the Spectra
The spectra is reduced with the Binospec pipeline described [here](https://bitbucket.org/chil_sai/binospec/wiki/Home). The 1D spectra is now split into individual fits file, one spectrum per file. They are distributed bundled up into tar files, whose contents are the 1D files with file names derived from object names submitted with the mask design. For the example here, obj_abs_1D.tar folder is the default folder used. The spectra from the pipeline still need further cleaning from the sky lines and bad data. This script produce new fits table with the cleaned spectrum for each file and preliminary plots for a quicklook.

```
python3 clean_spec.py
```

### 


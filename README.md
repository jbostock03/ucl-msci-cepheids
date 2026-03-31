# Calibrating the Cepheid Period-Luminosity Relation with _Gaia_ DR3
All code and files used in the Master's project for PHAS0097 Astrophysics Project at UCL.

## Contents
- [Code](https://github.com/jbostock03/ucl-msci-cepheids#Code)
- [Datasets](https://github.com/jbostock03/ucl-msci-cepheids#Datasets)
- [Extras](https://github.com/jbostock03/ucl-msci-cepheids#Extras)
- [Figures](https://github.com/jbostock03/ucl-msci-cepheids#Figures)
- [Acknowledgements](https://github.com/jbostock03/ucl-msci-cepheids#Acknowledgements)

## Code
Available are all .py and .pynb files used in the project:

**`selection_plotting.ipynb`**: _(main notebook)_ Full workflow of Bayesian distance estimation using a selection function _S_ = 1, then with MCMC sampling to deduce a sigmoid-type selection function, then redone with this. Also includes plotting for the report.

**`flat prior.ipynb`**: Same Notebook as above but with a flat selection function used in MCMC sampling: _S_ = _S0_

**`PL_relations.py`**:  Plots PL relation from least squares using either Gaia DR2 or DR3 data with 1/p as a distance estimator

**`plot_gaia_z.py`**: Plots the distribution of heights of Cepheids and a scatter plot of latitude v. longitude

**`plot_gaia_z_sns.py`**: Plots KDE of latitude v. longitude

**`riess_copy.py`**: Least squares fit using 1/p distance estimator and DR2 data

## Datasets
Available are all datasets used in the project in .csv format:

**`incl. l and b-result.csv`**: Full dataset of 13 077 DR3 Cepheids used including latitudes and longitudes

**`riess2018gaiadr2.csv`**: Dataset of 50 Cepheids as used by Riess et al. (2018b), including _HST_ photometry results and DR2 parallaxes

**`xmatch_riess_col_corr_gaia_source-result.csv`**: Cross-matched result from the _Gaia_ Archive between the above two datasets. ADQL queries used for these are in Appendix A in the report.

## Extras
Extra datasets used at points in the project are located in the folder **`/extra_data/`**:

**`init_cepheids.csv`**: Small sample of Cepheids not used in the final report

**`22jan_full_cepheids.csv`**: Incomplete sample of Cepheids not used in the final report

**`25feb_full_cepheids.csv`**: Incomplete sample of Cepheids not used in the final report

**`riess_col_corrected.csv`**: Formatted version of `riess2018gaiadr2.csv` for use in cross-matching

## Figures
All figures used in the report are located in the folder **`/figures/`**.

## Acknowledgements
Data: This work has made use of data from the European Space Agency (ESA) mission Gaia (https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing and Analysis Consortium (DPAC, https://www.cosmos.esa.int/web/gaia/dpac/consortium). Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.

Software: NumPy (Harris et al., 2020), pandas (McKinney et al., 2011), Matplotlib (Hunter, 2007), AstroPy (Astropy Collaboration et al., 2013), SciPy (Virtanen et al., 2020), seaborn (Waskom, 2021), emcee (Foreman-Mackey et al., 2013), corner.py (Foreman-Mackey et al., 2021), astroML (VanderPlas et al., 2012), gaiadr3 zeropoint (Lindegren et al., 2021).

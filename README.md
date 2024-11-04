# Affinity-Solvent-Effects
Used for analysis of boronic acid-diol titrations and their corresponding solvent-based affinity trends

## About
Author: Aidan Comiskey \
Contact: aidan.comiskey@columbia.edu

This repository contains raw UV/vis absorbance spectra, experimental parameters, and curvefitting parameters for a series of boronic acid-diol titrations to be published in a manuscript entitled "Solvent Effects in Boronic Acid-Diol Binding" by Aidan Comiskey (Columbia University) and Dr. Eric V Anslyn (Welch Regents Chair, Department of Chemistry, University of Texas at Austin) It also contains the code used to analyze this data, along with the output files, including plots and tables. The data analysis may be reproduced by running main.py.


## Data
The contents of the Data directory are as follows:
* Spectra - A directory containing unprocessed UV/vis spectra files (.csv) by titration
* direct_experimental_values.xlsx - Fitted values for direct titrations (i.e. between a boronic acid and indicator) to be used when fitting IDA titration isotherms
* experimental_parameters.csv - Information regarding how each titration experiment was conducted (concentrations used, solvent, identity of host and guest, etc.)
* fit_parameters.csv - Settings for fitting each isotherm (Whether indicator concentration was held constant, fitting algorithm, wavelengths analyzed, initial guess for Ka being fit and molar extinction coefficients, wavelengths analyzed, etc.)
* solvent_params.csv - Solvent parameters for analysis of solvent-mediated Ka trends

## Running the code
### Necessary Packages:
The following packages must be installed in the virtual environment for the code to be successfully executed: cycler, lmfit, matplotlib, numpy,  openpyxl, pandas, scipy, sklearn, statsmodels, and tqdm.

### Building the Necessary Directories:
To rerun the data analysis, the Data directory and all Python files must be downloaded. If output files (contained in the Plots, Saved Titration Fits, and Tables directories) are downloaded, they will be rewritten as the analysis is rerun. To rerun the analysis without downloading the output files, empty directories for the output files must be created within the working directory. The structure of these directories (and their subfolders) is given below, with a brief description of each element in parentheses:
* Plots
  * Combined (side-by-side titration spectra and isotherms)
  * Isotherms (titration isotherms)
  * Model (plots of multivariate regression of Ka against solvent parameters)
    * contour
* Spectra (titration spectra)
* Saved Titration Fits (pickled TitrationGroup() object - post fitting)
* Tables
  * Fit (fitted parameter values for direct and IDA titrations)
  * Model Summaries (full model statistics for the multivariate regressions of Ka against chosen solvent parameters)
  * Solvent Regression (summarized model statistics for the single and multivariate regressions of Ka against all solvent parameters, all combinations of two solvent parameters)

### Runtime Parameters
The data analysis can be rerun in two ways: with or without refitting all titration isotherms. To refit all isotherms, the REFIT runtime parameter in main.py must be set to True. To reanalyze the solvent parameter-affinity correlations and output the fitted parameter tables without refitting all titrations, REFIT should be set to False, and the SAVED_FITS_FILENAME parameter should be set to the saved TitrationGroup object containing existing fits (note that this will not reproduce the titration spectra and isotherm plots, which are created as fitting is performed). For the most fitting run corresponding to the analysis found in the manuscript, this is "Saved Titration Fits/2024-11-03 21:07:45.861605".

To rerun the analysis, main.py must be executed. Output files will then populate their corresponding directories (note that Plots will only be populated if titrations are refit). 


from automation import analyze_all_titrations, update_direct_exp_values_excel
import pickle

DIRECT_LIST = ['168', '169', '129', '174', '93', '139', '157', '159', '147', '149']
IDA_LIST = ['170', '171', '133', '135', '177a', '177b', '183a', '181a', '183b', '181b', '178a', '178b', '180a', '180b',
            '184a', '182a', '184b', '182b']
FULL_LIST = DIRECT_LIST + IDA_LIST

# Input files:
EXPERIMENTAL_PARAMS = 'Data/experimental_parameters.csv'
FIT_PARAMS = 'Data/fit_parameters.csv'
DIRECT_EXPERIMENTAL_VALUES = 'Data/direct_experimental_values.xlsx'
SPECTRUM_FORMAT = "Data/Spectra/AC-{}.csv"

# Runtime Parameters, set to preference:
PRINT_REPORT = True  # Prints fit reports as each titration is fitted
CORRECT_BASELINE = True  # Corrects titration spectrum baseline across all wavelengths such that Abs @ 750 nm = 0
REFIT = True  # Runs fitting process again; if False, previously-fitted values from SAVED_FITS FILENAME are used
SAVED_FITS_FILENAME = "Saved Titration Fits/2024-11-03 18:20:33.529232"

if REFIT:
    update_direct_exp_values_excel(direct_titrations=DIRECT_LIST, experimental_params_file=EXPERIMENTAL_PARAMS,
                                   fit_params_file=FIT_PARAMS, experimental_values_file=DIRECT_EXPERIMENTAL_VALUES,
                                   spectrum_file_format=SPECTRUM_FORMAT, print_report=PRINT_REPORT,
                                   correct_baseline=CORRECT_BASELINE)

    titrations = analyze_all_titrations(titration_list=FULL_LIST, experimental_params_file=EXPERIMENTAL_PARAMS,
                                        fit_params_file=FIT_PARAMS, experimental_values_file=DIRECT_EXPERIMENTAL_VALUES,
                                        spectrum_file_format=SPECTRUM_FORMAT, print_report=PRINT_REPORT,
                                        correct_baseline=CORRECT_BASELINE)
else:
    file = open(SAVED_FITS_FILENAME, "rb")
    titrations = pickle.load(file)
    file.close()

# Saves fit output table, solvent param regression value table, and full summary of multivariate for given predictors
titrations.get_fit_params_tables(direct_titration_list=DIRECT_LIST, IDA_titration_list=IDA_LIST)
titrations.get_regressions(standardize_multi=True)
titrations.get_model_summary(predictors=('Dielectric', 'beta'))

if REFIT:
    titrations.pickle()  # Saves fitted values (in titration object) for later analysis

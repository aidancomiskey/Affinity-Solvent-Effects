import pickle
import numpy as np
import pandas as pd

from lmfit import minimize, report_fit
from datetime import datetime
from fitter import construct_fit_params, global_residuals, model_absorbances_global
from plot import plot_spectra, plot_fit, plot_together, title_plot
from tables import get_fit_params_table, get_regressions, get_model_summary


# Custom object to load in a single titration experiment
# All dataframes (dfs) are loaded from TitrationGroup instance NOTE: TitrationGroup must be initialized first
class Titration:
    def __init__(self, experiment_number: str, experimental_params_df, fit_params_df, experimental_values_df,
                 spectrum_file_format):

        if not isinstance(experiment_number, str):
            raise TypeError("Experiment number must be in string format")

        # Initialize values to be fitted to None
        self.guest_concentrations_fit = None
        self.fit_output = None
        self.output_K = None
        self.fit_dict = None
        self.output_coefficients = None
        self.output_indicator_concentration = None

        # Load in parameters from experimental_params_df
        selected_experimental_params = experimental_params_df.loc[experiment_number]
        self.experiment_number = experiment_number
        self.solvent = selected_experimental_params["Solvent"]
        self.host = selected_experimental_params["host"]
        self.guest = selected_experimental_params["guest"]
        self.indicator = selected_experimental_params["indicator"]
        self.host_concentration = selected_experimental_params["host_concentration"]
        self.indicator_concentration = selected_experimental_params["indicator_concentration"]
        self.type = selected_experimental_params["type"]

        # Loading in fit parameters from fit_params_df
        selected_fit_params = fit_params_df.loc[experiment_number]
        self.lambda1 = selected_fit_params["lambda1"]
        self.lambda2 = selected_fit_params["lambda2"]
        self.indicator_flex = True if selected_fit_params["indicator_flex"] == "Flex" else False
        self.fit_method = selected_fit_params["fit_method"]
        self.k_guess = selected_fit_params["k_guess"]

        # For direct titrations: loads in fitting estimations for molar extinction coefficients
        if self.type == "Direct":
            self.guess_epsilon_host1 = selected_fit_params["guess_epsilon_host1"]
            self.guess_epsilon_complex1 = selected_fit_params["guess_epsilon_complex1"]
            self.guess_epsilon_host2 = selected_fit_params["guess_epsilon_host2"]
            self.guess_epsilon_complex2 = selected_fit_params["guess_epsilon_complex2"]
            self.experimental_k_i = None
        elif self.type == "IDA":
            associated_direct_titration = selected_fit_params["associated_direct_titration"]
            selected_experimental_values = experimental_values_df.loc[associated_direct_titration]
            print(selected_experimental_values)
            self.guess_epsilon_host1 = selected_experimental_values["epsilon_host1"]
            self.guess_epsilon_complex1 = selected_experimental_values["epsilon_complex1"]
            self.guess_epsilon_host2 = selected_experimental_values["epsilon_host2"]
            self.guess_epsilon_complex2 = selected_experimental_values["epsilon_complex2"]
            self.experimental_k_i = selected_experimental_values["experimental_k_i"]

        # Load in preformatted spectrum file and split into concentrations, wavelengths, and absorbances
        spectrum_file = spectrum_file_format.format(experiment_number)
        spectrum_data = np.genfromtxt(spectrum_file, skip_header=1, max_rows=602, delimiter=",", dtype=float)
        self.guest_concentrations = spectrum_data[0, 1:]
        self.wavelengths = spectrum_data[1:, 0]
        self.absorbances = spectrum_data[1:, 1:]

        # Select absorbances corresponding to selected wavelengths to construct experimental isotherms
        self.n_points = len(self.guest_concentrations)
        self.selected_absorbances = np.zeros((2, self.n_points))
        self.selected_absorbances[0] = self.absorbances[self.wavelengths == self.lambda1].flatten()
        self.selected_absorbances[1] = self.absorbances[self.wavelengths == self.lambda2].flatten()

    # Corrects baselines by shifting equally across spectrum, referenced to wavelength where abs should be 0
    def correct_baselines(self, ref_wavelength, ref_absorbance):
        correction_deltas = ref_absorbance - self.absorbances[self.wavelengths == ref_wavelength].flatten()
        self.absorbances += correction_deltas
        self.selected_absorbances[0] = self.absorbances[self.wavelengths == self.lambda1].flatten()
        self.selected_absorbances[1] = self.absorbances[self.wavelengths == self.lambda2].flatten()

    # Creates plot of spectra over the course of the titration
    def plot_spectra(self):
        plot_spectra(self)

    # Fits isotherm curves at two wavelengths using global analysis, nonlinear methods
    def fit(self, print_report=True):
        fit_params = construct_fit_params(self)
        self.fit_output = minimize(global_residuals, params=fit_params, args=(self.guest_concentrations, self),
                                   method=self.fit_method)
        if print_report:
            print(title_plot(self))
            report_fit(self.fit_output)

        # Records output_K (affinity as found by fitter)
        if self.type == "Direct":
            self.output_K = self.fit_output.params["k_i"].value
        elif self.type == "IDA":
            self.output_K = self.fit_output.params["k_g"].value

        # Records indicator concentration after fitting (=/= input concentration if this value is fit)
        self.output_indicator_concentration = self.fit_output.params['indicator_concentration'].value

        # Records output molar extinction coefficients at wavelengths 1 and 2
        self.output_coefficients = {"epsilon_host1": self.fit_output.params["epsilon_host1"].value,
                                    "epsilon_host2": self.fit_output.params["epsilon_host2"].value,
                                    "epsilon_complex1": self.fit_output.params["epsilon_complex1"].value,
                                    "epsilon_complex2": self.fit_output.params["epsilon_complex2"].value}

        # Gets range of concentrations for continuous plot of isotherm fit
        fit_concentrations = np.linspace(start=self.guest_concentrations[0],
                                         stop=self.guest_concentrations[-1],
                                         num=100)
        # Models absorbance values for isotherm fit
        fit_absorbances = model_absorbances_global(params=self.fit_output.params,
                                                   guest_concentrations=fit_concentrations,
                                                   titration=self)

        self.fit_dict = {"concentrations": fit_concentrations,
                         "absorbances": fit_absorbances}  # Packing conc. range and fit abs into dict

    # Plots isotherm with line of best fit
    def plot_fit(self):
        plot_fit(self)

    # Plots isotherm alongside full spectra over the course of the titration
    def plot_together(self):
        plot_together(self)

    # Packaging fitting and plotting functions into single analysis function
    def analyze_titration(self, print_report=True):
        self.plot_spectra()
        self.fit(print_report=print_report)
        self.plot_fit()
        self.plot_together()


# Custom object for group of titrations (i.e. titrations completed in different solvents to analyze solvent trends)
class TitrationGroup:

    # Initializes Titration Group with files describing experimental parameters, guess values for Ka and molar
    # extinction coefficients, and the directory path to experimental spectra
    def __init__(self, experimental_params_file, fit_params_file, experimental_values_file, spectrum_file_format):
        # Loads experimental parameters
        self.experimental_params_df = pd.read_csv(experimental_params_file, index_col=0)
        self.experimental_params_df.index = self.experimental_params_df.index.astype(str)

        # Loads fitting parameters
        self.fit_params_df = pd.read_csv(fit_params_file, index_col=0)
        self.fit_params_df.index = self.fit_params_df.index.astype(str)

        self.spectrum_file_format = spectrum_file_format  # Formatting for how spectra files are named

        # Loads in experimental values from direct titrations
        self.experimental_values_filename = experimental_values_file
        self.experimental_values_df = pd.read_excel(self.experimental_values_filename, index_col=0)

        # Initialize empty/0 attributes that will be filled upon completing fits for all titrations
        self.affinities = {}  # Stores affinities for each experiment number in a dictionary
        self.titrations = {}

    def __getitem__(self, item):
        return self.titrations[item]

    # Initialize Titration object as an attribute of the TitrationGroup.titrations dict
    def load_titration(self, experiment_number):
        self.titrations[experiment_number] = Titration(experiment_number=experiment_number,
                                                       experimental_params_df=self.experimental_params_df,
                                                       fit_params_df=self.fit_params_df,
                                                       experimental_values_df=self.experimental_values_df,
                                                       spectrum_file_format=self.spectrum_file_format)

    # Organizes output affinities (i.e. fit Ka values) for all Titrations into dictionary
    def get_affinities(self):
        for titration in self.titrations:
            self.affinities[titration] = self[titration].output_K

    # Store titration group after fitting (to work with results without repeating fitting routine)
    def pickle(self):
        current_time = datetime.now()
        filename = f'Saved Titration Fits/{current_time}'
        file = open(filename, 'wb')
        pickle.dump(obj=self, file=file)
        file.close()
        print(f'TitrationGroup saved!          See {filename}')

    # Save fitted parameters to table
    def get_fit_params_tables(self, direct_titration_list, IDA_titration_list):
        get_fit_params_table(titration_group=self, experiment_list=direct_titration_list, IDA=False)
        get_fit_params_table(titration_group=self, experiment_list=IDA_titration_list, IDA=True)

    # Run statistical analysis on titration group
    def get_regressions(self, standardize_multi=True):
        get_regressions(titration_group=self, standardize_multi=standardize_multi)

    # Get model summary from predictors
    def get_model_summary(self, predictors):
        get_model_summary(titration_group=self, predictors=predictors, plot=True)

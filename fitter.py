import numpy as np
from lmfit import Parameters
from uvvis_models import model_absorbances_direct, model_absorbances_IDA

K_I_FLEX = 0.2
COEFF_FLEX = 0.4
STRICT_COEFF_BOUNDS = False
VARY_HOST_CONC = False
CONC_FLEX = 0.1


# Determines fitting bounds for molar extinction coefficients
def get_coeff_bounds(value):
    if STRICT_COEFF_BOUNDS:
        max = (1 + COEFF_FLEX) * value if value > 80 else 112
        min = (1 - COEFF_FLEX) * value if value > 100 else 0
    else:
        max = 20000
        min = 0
    bounds = {'max': max, 'min': min}
    return bounds


# Returns lmfit.Parameters() object for a given titration
def construct_fit_params(titration):
    fit_params = Parameters()
    fit_params.add('indicator_concentration', value=titration.indicator_concentration,
                   vary=titration.indicator_flex,
                   min=(1 - CONC_FLEX) * titration.indicator_concentration,
                   max=(1 + CONC_FLEX) * titration.indicator_concentration)
    if titration.type == "Direct":
        fit_params.add('epsilon_host1', value=titration.guess_epsilon_host1, min=0, max=20000)
        fit_params.add('epsilon_host2', value=titration.guess_epsilon_host2, min=0, max=20000)
        fit_params.add('epsilon_complex1', value=titration.guess_epsilon_complex1, min=0, max=20000)
        fit_params.add('epsilon_complex2', value=titration.guess_epsilon_complex2, min=0, max=20000)
        fit_params.add('k_i', value=titration.k_guess, min=0, max=1e6)
    elif titration.type == "IDA":
        fit_params.add('epsilon_host1', value=titration.guess_epsilon_host1,
                       min=get_coeff_bounds(titration.guess_epsilon_host1)['min'],
                       max=get_coeff_bounds(titration.guess_epsilon_host1)['max'])
        fit_params.add('epsilon_host2', value=titration.guess_epsilon_host2,
                       min=get_coeff_bounds(titration.guess_epsilon_host2)['min'],
                       max=get_coeff_bounds(titration.guess_epsilon_host2)['max'])
        fit_params.add('epsilon_complex1', value=titration.guess_epsilon_complex1,
                       min=get_coeff_bounds(titration.guess_epsilon_complex1)['min'],
                       max=get_coeff_bounds(titration.guess_epsilon_complex1)['max'])
        fit_params.add('epsilon_complex2', value=titration.guess_epsilon_complex2,
                       min=get_coeff_bounds(titration.guess_epsilon_complex2)['min'],
                       max=get_coeff_bounds(titration.guess_epsilon_complex2)['max'])
        fit_params.add('k_i', value=titration.experimental_k_i,
                       min=(1 - K_I_FLEX) * titration.experimental_k_i,
                       max=(1 + K_I_FLEX) * titration.experimental_k_i)
        fit_params.add('host_concentration', value=titration.host_concentration, vary=VARY_HOST_CONC,
                       min=(1 - CONC_FLEX) * titration.host_concentration,
                       max=(1 + CONC_FLEX) * titration.host_concentration)
        fit_params.add('k_g', value=titration.k_guess, min=0, max=1e6)
    return fit_params


# Returns modeled absorbances at both selected wavelengths given inputted parameter values
def model_absorbances_global(params, guest_concentrations, titration):
    model_absorbances = []
    if titration.type == "Direct":
        for i, wavelength in enumerate([titration.lambda1, titration.lambda2], start=1):
            model_absorbances.append(model_absorbances_direct(guest_concentrations=guest_concentrations,
                                                              indicator_concentration=params["indicator_concentration"],
                                                              k_i=params["k_i"],
                                                              epsilon_host=params[f"epsilon_host{i}"],
                                                              epsilon_complex=params[f"epsilon_complex{i}"]))
    elif titration.type == "IDA":
        for i, wavelength in enumerate([titration.lambda1, titration.lambda2], start=1):
            model_absorbances.append(model_absorbances_IDA(guest_concentrations=guest_concentrations,
                                                           indicator_concentration=params["indicator_concentration"],
                                                           k_i=params["k_i"],
                                                           epsilon_host=params[f"epsilon_host{i}"],
                                                           epsilon_complex=params[f"epsilon_complex{i}"],
                                                           host_concentration=params["host_concentration"],
                                                           k_g=params["k_g"]))
    return np.array(model_absorbances)


# Returns flattened array of residuals for modeled absorbance values at two wavelengths
def global_residuals(fit_params, guest_concentrations, titration):
    selected_absorbances_flattened = titration.selected_absorbances.flatten()
    model_absorbances = model_absorbances_global(params=fit_params,
                                                 guest_concentrations=guest_concentrations,
                                                 titration=titration)
    model_absorbances_flattened = model_absorbances.flatten()
    residuals_flattened = selected_absorbances_flattened - model_absorbances_flattened
    return residuals_flattened

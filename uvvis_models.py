import numpy as np


# Models absorbance values for 1:1 UV/vis isotherm, equation adapted from Hargrove et al.
def model_absorbances_direct(guest_concentrations, indicator_concentration, k_i, epsilon_host, epsilon_complex):
    x = guest_concentrations
    i0 = indicator_concentration
    complex_concentrations = 0.5 * ((i0 + (1 / k_i) + x) - np.sqrt((i0 + (1 / k_i) + x) ** 2 - 4 * i0 * x))
    model_absorbances = epsilon_host * (i0 - complex_concentrations) + epsilon_complex * complex_concentrations
    return model_absorbances


# Models absorbance values for indicator displacement isotherm, equation adapted from Hargrove et al.
def model_absorbances_IDA(guest_concentrations, indicator_concentration, k_i, epsilon_host, epsilon_complex,
                          host_concentration, k_g):
    x = guest_concentrations
    h0 = host_concentration
    i0 = indicator_concentration
    a = k_i * k_g
    b = k_i + k_g + k_i * k_g * i0 + k_i * k_g * x - k_i * k_g * h0
    c = 1 + k_i * i0 + k_g * x - k_i * h0 - k_g * h0
    d = -h0
    free_host_concentrations = np.array([np.polynomial.polynomial.Polynomial([d, y, x, a]).roots()
                                         for x, y in zip(b, c)])[:, 2]
    model_absorbances = i0 * (epsilon_host + epsilon_complex * k_i * free_host_concentrations) / (
            1 + k_i * free_host_concentrations)
    return model_absorbances

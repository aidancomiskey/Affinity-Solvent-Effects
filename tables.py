import pandas as pd
from solvent_correlation import single_regression, multi_regression

SINGLE_REGRESSORS = ["Dielectric", "Z", "ET30", "Pi*", "alpha", "beta"]

# Combines SINGLE_REGRESSORS into pairs
REGRESSOR_PAIRS = []
for i, r1 in enumerate(SINGLE_REGRESSORS[:-1], 1):
    for r2 in SINGLE_REGRESSORS[i:]:
        REGRESSOR_PAIRS.append((r1, r2))


# Returns table of all fitted parameters by titration experiment
def get_fit_params_table(experiment_list, titration_group, IDA=False):
    output_dict = {"exp_number": [], "free_epsilons": [], "complex_epsilons": [], "Ks_indicator": [], "i_conc_fit": []}
    if IDA:
        output_dict["Ks_guest"] = []

    for experiment_number in experiment_list:
        titration = titration_group[experiment_number]
        output_params = titration.fit_output.params
        output_dict["exp_number"].extend((experiment_number, ""))
        output_dict["free_epsilons"].extend(
            (output_params["epsilon_host1"].value, output_params["epsilon_host2"].value))
        output_dict["complex_epsilons"].extend(
            (output_params["epsilon_complex1"].value, output_params["epsilon_complex2"].value))
        output_dict["Ks_indicator"].extend((output_params["k_i"].value, ""))
        if IDA:
            output_dict["Ks_guest"].extend((output_params["k_g"].value, ""))
        output_dict["i_conc_fit"].extend((output_params["indicator_concentration"].value, ""))

    df_output = pd.DataFrame(output_dict)
    name = "IDA" if IDA else "Direct"
    fit_summary_output_filename = f"Tables/Fit/{name}.xlsx"
    df_output.to_excel(fit_summary_output_filename, index=False)
    print("See fit summaries at ", fit_summary_output_filename)


# Returns table detailing single solvent parameter-affinity regression results for all predictors by boronic
# acid-diol pair
def single_regression_table(titration_group, print_results=False):
    dfs = [single_regression(titration_group=titration_group, predictor=predictor, print_results=print_results) for
           predictor in SINGLE_REGRESSORS]
    results_df = pd.concat(dfs)
    return results_df


# Returns table detailing multivariate solvent parameter-affinity regression results for all predictor combinations
# by boronic acid-diol pair
def multivariate_regression_table(titration_group, print_results=False, standardize=True):
    dfs = [multi_regression(titration_group=titration_group, standardize=standardize, predictors=predictors,
                            print_results=print_results) for predictors in REGRESSOR_PAIRS]
    results_df = pd.concat(dfs)
    return results_df


def model_summary_table(predictors, titration_group, print_results=False):
    df = multi_regression(titration_group=titration_group,
                          predictors=predictors,
                          report_full_model=True,
                          print_results=print_results,
                          standardize=True)
    return df


def get_regressions(titration_group, standardize_multi=True, tables_filepath="Tables/Solvent Regression/"):
    single_regression_df = single_regression_table(titration_group=titration_group)
    single_regression_filename = tables_filepath + "Single.xlsx"
    single_regression_df.to_excel(single_regression_filename)
    print("See Single Variable Regression Summaries at ", single_regression_filename)

    multi_regression_df = multivariate_regression_table(titration_group=titration_group, standardize=standardize_multi)
    multi_regression_filename = tables_filepath + "Multivariate.xlsx"
    multi_regression_df.to_excel(multi_regression_filename)
    print("See Multivariable Regression Summaries at ", multi_regression_filename)


# Get summary table and contour plots for model by predictors:
def get_model_summary(titration_group, predictors, plot=False):
    model_summary_df = model_summary_table(predictors=predictors, titration_group=titration_group)
    model_summary_filename = f"Tables/Model Summaries/{predictors[0]} & {predictors[1]}.xlsx"
    model_summary_df.to_excel(model_summary_filename)
    print("See model summary at ", model_summary_filename)

    # Get 3D plots of model
    if plot:
        for style in ["contour"]:
            multi_regression(titration_group=titration_group, predictors=predictors,
                             plot_3D=True, plot_style=style, standardize=False)

import pandas as pd
import statsmodels.api as ssm
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor
from plot import plot_multivariate
from sklearn.preprocessing import StandardScaler

SOLVENT_PARAMETERS = "Data/solvent_params.csv"
BORONIC_ACIDS = ["PBA", "ompPBA"]
DIOLS = ["Alizarin", "Catechol", "Hydrobenzoin"]
BINDING_PAIRS = [(boronic, diol) for boronic in BORONIC_ACIDS for diol in DIOLS]
EXPERIMENT_LISTS = [["129", "93", "157", "147"], ["170", "133", "183a", "178a", "184a"],
                    ["135", "181a", "178b", "182a"], ["174", "139", "159", "149"],
                    ["171", "177a", "183b", "180a", "184b"], ["177b", "181b", "180b", "182b"]]

solvent_param_df = pd.read_csv(SOLVENT_PARAMETERS, index_col=0)


# Performs linear regressions by boronic acid-diol pair for a single solvent parameter (predictor)
def single_regression(titration_group, predictor, print_results=False):
    result_dict = {}
    for pair, experiments in zip(BINDING_PAIRS, EXPERIMENT_LISTS):
        titrations = [titration_group[exp] for exp in experiments]
        solvents = [titration.solvent for titration in titrations]
        K_values = [titration.output_K for titration in titrations]
        regression_df = pd.DataFrame({"solvent": solvents, "K_values": K_values})
        regression_df[predictor] = regression_df["solvent"].map(solvent_param_df[predictor])
        result = stats.linregress(x=regression_df[predictor], y=regression_df["K_values"])
        result_dict[pair] = {f"coefficient_{predictor}": result.slope,
                             f"r_squared_{predictor}": result.rvalue ** 2,
                             f"p_value_{predictor}": result.pvalue}  # Choose which info to report in Dataframe
        if print_results:
            print(pair)
            print(pd.DataFrame(result_dict))
    return pd.DataFrame(result_dict)


# Performs multivariate linear regressions by boronic acid-diol pair for two solvent parameters (predictors)
def multi_regression(titration_group, predictors, standardize=True, report_full_model=False, plot_3D=False,
                     plot_style="3d", print_results=False):
    p1, p2 = predictors
    result_dict = {}
    for pair, experiments in zip(BINDING_PAIRS, EXPERIMENT_LISTS):
        titrations = [titration_group[exp] for exp in experiments]
        solvents = [titration.solvent for titration in titrations]
        K_values = [titration.output_K for titration in titrations]
        regression_df = pd.DataFrame({"solvent": solvents, "K_values": K_values})
        regression_df[p1] = regression_df["solvent"].map(solvent_param_df[p1])
        regression_df[p2] = regression_df["solvent"].map(solvent_param_df[p2])

        y = regression_df["K_values"]
        X = regression_df[[p1, p2]]

        # Standardize X data by converting to Z score
        if standardize:
            scaler = StandardScaler()
            X_transform = pd.DataFrame(columns=[p1, p2], index=X.index)
            X_transform[p1] = scaler.fit_transform(X[p1].to_numpy().reshape(-1, 1))
            X_transform[p2] = scaler.fit_transform(X[p2].to_numpy().reshape(-1, 1))
            X = ssm.add_constant(X_transform)

        else:
            X = ssm.add_constant(X)

        result = ssm.OLS(y, X).fit()
        vif = variance_inflation_factor(X, 1)
        if report_full_model:
            result_dict[pair] = {f"coefficient ({p1})": result.params[p1],
                                 f"Std. Error ({p1})": result.HC1_se[p1],
                                 f"p-value ({p1})": result.pvalues[p1],
                                 f"coefficient ({p2})": result.params[p2],
                                 f"Std. Error ({p2})": result.HC1_se[p2],
                                 f"p-value ({p2})": result.pvalues[p2],
                                 "Intercept Value": result.params["const"],
                                 "Intercept Std. Error": result.HC1_se["const"],
                                 "Intercept p-value": result.pvalues["const"],
                                 "n =": result.nobs,
                                 "VIF": vif,
                                 "R-squared": result.rsquared,
                                 "R-squared_adj": result.rsquared_adj,
                                 "F-statistic": result.fvalue,
                                 "p(F-statistic)": result.f_pvalue,
                                 "Log-Likelihood": result.llf,
                                 "AIC": result.aic,
                                 "BIC": result.bic}  # Choose which info to report in Dataframe
        else:
            result_dict[pair] = {f"r_squared_{p1}_{p2}": result.rsquared,
                                 f"p-value ({p1})": result.pvalues[p1],
                                 f"p-value ({p2})": result.pvalues[p2],
                                 f"VIF": vif}  # Choose which info to report in Dataframe
        if plot_3D:
            reduced_y = y/1e4
            title = f"{pair[0]} & {pair[1]}"
            plot_multivariate(X=X,
                              y=reduced_y,
                              p1=p1,
                              p1_coeff=result.params[p1] / 1e4,
                              p2=p2,
                              p2_coeff=result.params[p2] / 1e4,
                              intercept=result.params["const"] / 1e4,
                              y_label="Affinity (×10⁴)",
                              title=title,
                              filename=f"Plots/Model/{plot_style}/{title}.png",
                              style=plot_style,
                              standardize=standardize)
        if print_results:
            print(f"Pair: {pair[0]} & {pair[1]}")
            print(result.summary())
    return pd.DataFrame(result_dict)

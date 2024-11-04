import numpy as np
from cycler import cycler
import matplotlib

matplotlib.use('TkAgg')

from matplotlib import pyplot as plt

SUPERSCRIPT_MAP = {"0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹"}

# Set for spectra gradient; format (r,g,b) where each value is between 0 and 1
color_start = (0.88, 0.12, 0.12)  # red
color_end = (0.12, 0.12, 0.88)  # blue


def title_plot(titration):
    title = None
    if titration.type == "Direct":
        title = f'Direct Titration of {titration.indicator} & {titration.guest} in {titration.solvent} \n'
    elif titration.type == "IDA":
        title = f'IDA Titration of {titration.indicator}:{titration.host} complex & {titration.guest} in ' \
                f'{titration.solvent} \n'
    return title


def annotate_K(titration):
    subscript = "I" if titration.type == "Direct" else "G"
    value = f'{titration.output_K:.2E}'
    sigfigs, exponent = value.split("E+0")
    text_K = '$K_{{{subscript}}} = {sigfigs}×10{exponent} M^{{-1}}$'.format(subscript=subscript,
                                                                            sigfigs=sigfigs,
                                                                            exponent=SUPERSCRIPT_MAP[exponent])
    return text_K


def annotate_conc(titration):
    text = f'[{titration.indicator}] = {titration.indicator_concentration * 1000000:.0f} $\mu$M'
    if titration.type == "IDA":
        text += f' \n [{titration.host}] = {titration.host_concentration * 1000000:.0f} $\mu$M'
    return text


# Returns color gradient for plotting titration spectra
def get_gradient(start, end, steps):
    grad = np.linspace(0, 1, steps)
    color_diff = np.array(end) - np.array(start)
    diffs = np.outer(grad, color_diff)
    colors = diffs + np.array(start)
    return list(map(tuple, colors))


# Plots all spectra for a single titration
def plot_spectra(titration):
    custom_cycler_spectra = cycler(color=get_gradient(color_start, color_end, titration.n_points))
    fig, ax = plt.subplots()
    ax.set_prop_cycle(custom_cycler_spectra)
    ax.plot(titration.wavelengths, titration.absorbances, linewidth=0.5)
    ax.set(xlim=[200, 900], ylim=[0, 1.2], xlabel='Wavelength (nm)', ylabel='Absorbance', title=title_plot(titration))
    plt.savefig(f'Plots/Spectra/{titration.experiment_number}.svg')
    plt.close()


# Plots line of best fit and experimental isotherm for a given titration
def plot_fit(titration):
    custom_cycler_isotherm = cycler(color=["k", "r", "k", "r"])
    fig, ax = plt.subplots()
    ax.set_prop_cycle(custom_cycler_isotherm)
    x = titration.guest_concentrations
    y = titration.selected_absorbances.T
    ax.plot(x, y, ".", label=[f'{int(titration.lambda1)} nm', f'{int(titration.lambda2)} nm'])
    ax.plot(titration.fit_dict["concentrations"], titration.fit_dict["absorbances"].T)
    ax.set(xlabel=f'{titration.guest} Concentration [M]', ylabel='Absorbance', title=title_plot(titration))
    bottom, top = ax.get_ylim()
    ax.set_ylim(top=1.3 * top)
    ax.legend(loc="upper right")
    ax.text(0.8, 0.8, annotate_K(titration), horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    plt.savefig(f'Plots/Isotherms/{titration.experiment_number}.svg')
    plt.close()


# Plots titration spectra and isotherm side by side
def plot_together(titration):
    # Getting cyclers for spectra and isotherm colors
    custom_cycler_spectra = cycler(color=get_gradient(color_start, color_end, titration.n_points))
    custom_cycler_isotherm = cycler(color=["k", "r", "k", "r"])
    # Create plot
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=False, figsize=(8.5, 4))
    # Set colors and plot spectra
    ax1.set_prop_cycle(custom_cycler_spectra)
    ax1.plot(titration.wavelengths, titration.absorbances, linewidth=0.5)
    ax1.set(xlim=[200, 900], ylim=[0, 1.2], xlabel='Wavelength (nm)', ylabel='Absorbance')
    ax1.text(0.95, 0.95, annotate_conc(titration), horizontalalignment="right", verticalalignment="top",
             transform=ax1.transAxes)

    # Plot Isotherm:
    x = titration.guest_concentrations * 1000  # Convert guest concentrations to mM values
    y = titration.selected_absorbances.T  # Transpose isotherm values
    ax2.set_prop_cycle(custom_cycler_isotherm)
    ax2.plot(x, y, ".", label=[f'{int(titration.lambda1)} nm', f'{int(titration.lambda2)} nm'])  # Plot points
    ax2.plot(titration.fit_dict["concentrations"] * 1000, titration.fit_dict["absorbances"].T,
             linewidth=0.7)  # Plot fit
    ax2.set(xlabel=f'{titration.guest} Concentration [mM]', ylabel='Absorbance')

    bottom, top = ax2.get_ylim()  # only keep if sharey=False
    ax2.set_ylim(top=1.3 * top)  # only keep if sharey=False

    ax2.legend(loc="upper right", frameon=False)
    ax2.text(0.05, 0.95, annotate_K(titration), horizontalalignment='left', verticalalignment='top',
             transform=ax2.transAxes)

    # Title subplots:
    ax1.set_title("(a) Spectra", loc="left")
    ax2.set_title("(b) Isotherm", loc="left")

    fig.tight_layout(rect=(0, 0, 1, 1))

    plt.savefig(f'Plots/Combined/{titration.experiment_number}.png', dpi=300)
    plt.close()


# Plots multivariate model relating solvent parameters to affinity values for a given boronic acid-diol pair
def plot_multivariate(X, y, p1, p1_coeff, p2, p2_coeff, intercept, y_label, title, filename, style, standardize):
    if standardize:
        x1dim = np.linspace(-2, 2, 100)
        x2dim = np.linspace(-1.5, 1.5, 100)
    else:
        x1dim = np.linspace(20, 80, 100)
        x2dim = np.linspace(0.2, 0.8, 100)
    meshx1, meshx2 = np.meshgrid(x1dim, x2dim)
    predict_y = (p1_coeff * meshx1 + p2_coeff * meshx2 + intercept)

    if style == "3d":
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.set_title(title)
        ax.scatter(X[p1], X[p2], y)
        ax.plot_surface(meshx1, meshx2, predict_y, alpha=0.2, color=[0, 1, 0])
        ax.set_xlim(20, 80)
        ax.set_ylim(0.2, 0.8)
        ax.set_zlim(0, 25)
        ax.set_xlabel(p1)
        ax.set_ylabel(p2)
        ax.set_zlabel(y_label)
        plt.show()

    elif style == "contour":
        fig = plt.figure(figsize=(4, 3))
        ax = fig.add_subplot()

        contour = ax.contourf(meshx1, meshx2, predict_y, levels=50)
        fig.colorbar(contour, ax=ax, label=y_label)

        ax.scatter(X[p1], X[p2], c=y, s=80, edgecolors="k", norm=contour.norm)

        ax.set_xlabel(p1, labelpad=10)
        p2_label = "β" if p2 == "beta" else p2
        ax.set_ylabel(p2_label, rotation="horizontal", labelpad=15)
        ax.set_title(title, pad=10)

        fig.tight_layout(rect=(0, 0, 1, 1))

        plt.savefig(filename, dpi=300)
        plt.close()

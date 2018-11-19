from __future__ import absolute_import, division, print_function, unicode_literals

import uncertainpy as un
import chaospy as cp
import numpy as np
import matplotlib.pyplot as plt
import string
import os

from burstiness import is_bursting, is_regular, is_spiking, min_spike_amplitude

from uncertainpy.plotting.prettyplot.prettyplot import prettyBar, set_latex_font, set_style, spines_color, prettyPlot


# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
figure_format = ".eps"
figure_folder = "figures"


# Simulation parameters
discard = 10000                # ms
simulation_time = 60000        # ms
noise_amplitude = 0            # mV
stimulus_amplitude = 0         # nA

# Uncertainty quantification parameters
polynomial_order = 5
features_to_run = ["is_regular", "is_bursting", "is_spiking"]



def plot_sobol_feature(data, feature, ax):
    """
    Plot the total-order Sobol indices of `feature` from `data` on `ax`.

    Parameters
    ----------
    data : uncertainpy.Data
        A data object that contains the results from the uncertainty quantification.
    feature : str
        Name of the feature to plot.
    ax : matplotlib ax Object
        The ac Object to plot on.
    """
    style = "seaborn-darkgrid"

    width = 0.2
    index = np.arange(1, len(data.uncertain_parameters)+1)*width


    latex_labels = {"g_K": r"$g_\mathrm{K}$",
                    "g_Ca": r"$g_\mathrm{Ca}$",
                    "g_SK": r"$g_\mathrm{SK}$",
                    "g_Na": r"$g_\mathrm{Na}$",
                    "g_l": r"$g_\mathrm{l}$",
                    "g_BK": r"$g_\mathrm{BK}$"
    }

    xlabels = []
    for label in data.uncertain_parameters:
        xlabels.append(latex_labels[label])

        if feature.lower()  in ["rat", "medaka_1", "medaka_2"]:
            title = "Membrane potential"
        else:
            title = feature.replace("_", " ")
            title = title[0].upper() + title[1:]

        sensitivity = data[feature].sobol_total_average

        mean = data[feature].mean
        std = np.sqrt(data[feature].variance)
        unit = data[feature].labels

        if std < 1e-14:
            sensitivity = [0]*len(data.uncertain_parameters)

        if len(unit) > 0 and "(" in unit[0]:
            unit = unit[0].split("(")[-1].strip(")")
        else:
            unit = ""

        if feature == "spike_rate":
            mean *= 1000
            std *= 1000
            unit = "Hz"

        prettyBar(sensitivity,
                  xlabels=xlabels,
                  nr_colors=len(data.uncertain_parameters),
                  index=index,
                  ax=ax,
                  style=style)

        for tick in ax.get_xticklabels():
            tick.set_rotation(-40)

        ax.set_ylim([0, 1.15])
        ax.set_title(title, fontsize=titlesize)
        ax.text(0.25, 0.9,
                "Mean = {mean:.2{c}} {unit}".format(mean=mean,
                                                    c="e" if (abs(mean) < 1e-2 and mean != 0) else "f",
                                                    unit=unit),
                transform=ax.transAxes, fontsize=fontsize)
        ax.text(0.25, 0.78,
                "Std. = {std:.2{c}} {unit}".format(std=std,
                                                   c="e" if (abs(mean) < 1e-2 and mean != 0) else "f",
                                                   unit=unit),
                transform=ax.transAxes, fontsize=fontsize)

        ax.tick_params(labelsize=fontsize)




def plot_compare(rat, medaka_1, medaka_2):
    """
    Plot the total-order Sobol indices for the three models, rat,
    medaka 1 and medaka 2.

    Parameters
    ----------
    rat : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the rat model.
    medaka_1 : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the medaka 1 model.
    medaka_2 : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the medaka 2 model.
    """
    style = "seaborn-darkgrid"
    set_style(style)
    set_latex_font()

    plt.rcParams.update({"axes.titlepad": 12})
    fig, axes = plt.subplots(nrows=3, ncols=3, squeeze=False, sharey="row",
                             figsize=(figure_width, figure_width*0.8))



    set_style("seaborn-white")
    ax = fig.add_subplot(111, zorder=-10)
    spines_color(ax, edges={"top": "None", "bottom": "None",
                            "right": "None", "left": "None"})
    ax.tick_params(top=False, bottom=False, left=False, right=False, labelsize=fontsize,
                   labelbottom=False, labelleft=False)
    ax.set_ylabel('Total-order Sobol indices', labelpad=65, fontsize=labelsize)


    features = ["is_bursting", "is_regular", "is_spiking"]
    datas = [rat, medaka_1, medaka_2]

    yticks = np.arange(0, 1.1, 0.25)

    for i in range(3):
        for j in range(3):
            ax = axes[i][j]

            plot_sobol_feature(datas[i], features[j], ax)

            ax.text(-0.17, 1.08,
                    "\\textbf{{{}}}".format(string.ascii_uppercase[i] + str(j + 1)),
                    transform=ax.transAxes,
                    fontsize=titlesize)
            ax.set_title("")
            ax.set_yticks(yticks)


    axes[0][0].set_title("Isbursting")
    axes[0][1].set_title("Isregular")
    axes[0][2].set_title("Isfiring")

    axes[0][0].set_ylabel("RAT", labelpad=12)
    axes[1][0].set_ylabel("MEDAKA 1", labelpad=12)
    axes[2][0].set_ylabel("MEDAKA 2", labelpad=12)

    plt.tight_layout()
    plt.subplots_adjust(left=0.17)
    plt.savefig(os.path.join(figure_folder, "sensitivity" + figure_format))

    plt.rcdefaults()



def uq_rat():
    """
    Perform the uncertainty quantification and sensitivity analysis of the
    rat model.

    Returns
    -------
    data : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the rat model.
    """
    parameters = {"g_K": 9.55e-4,
                  "g_Ca": 6.34e-4,
                  "g_SK": 6.34e-4,
                  "g_l": 6.34e-5,
                  "g_BK": 2.4e-4} # Temporary value only

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))

    # Set full g_BK distribution
    parameters["g_BK"].distribution = cp.Uniform(0, 3.2e-4)

    # Initialize the features
    features = un.SpikingFeatures(new_features=[is_bursting, is_regular, is_spiking],
                                  features_to_run=features_to_run,
                                  logger_level="error",
                                  strict=False,
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False,
                                  min_amplitude=min_spike_amplitude)


    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="rat",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       plot=None,
                       filename="rat",
                       polynomial_order=polynomial_order,
                       save=False)

    return data


def uq_medaka_1():
    """
    Perform the uncertainty quantification and sensitivity analysis of the
    medaka 1 model.

    Returns
    -------
    data : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the medaka 1 model.
    """
    parameters = {"g_K": 9.55e-4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": 2.4e-4}  # Temporary value only

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))

    # Set full g_BK distribution
    parameters["g_BK"].distribution = cp.Uniform(0, 3.2e-4)

    # Initialize the features
    features = un.SpikingFeatures(new_features=[is_bursting, is_regular, is_spiking],
                                  features_to_run=features_to_run,
                                  logger_level="error",
                                  strict=False,
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False,
                                  min_amplitude=min_spike_amplitude)

    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka_1",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       plot=None,
                       filename="medaka_1",
                       polynomial_order=polynomial_order,
                       save=False)

    return data






def uq_medaka_2():
    """
    Perform the uncertainty quantification and sensitivity analysis of the
    medaka 2 model.

    Returns
    -------
    data : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the medaka 2 model.
    """
    parameters = {"g_K": 9.55e-4*1.4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4*3,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": 3.2e-4*4*0.67}  # Temporary value only

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))

    # Set full g_BK distribution
    parameters["g_BK"].distribution = cp.Uniform(0, 3.2e-4*4)

    # Initialize the features
    features = un.SpikingFeatures(new_features=[is_bursting, is_regular, is_spiking],
                                  features_to_run=features_to_run,
                                  strict=False,
                                  logger_level="error",
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False,
                                  min_amplitude=min_spike_amplitude)

    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka_2",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       filename="medaka_2",
                       plot=None,
                       polynomial_order=polynomial_order,
                       save=False)

    return data


def compare():
    """
    Perform the uncertainty quantification and sensitivity analysis of the
    rat, medaka 1 and medaka 2 models and compare the total-order Sobol
    indices of the .

    Returns
    -------
    data : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the rat model.
    """
    data_rat = uq_rat()
    data_medaka_1 = uq_medaka_1()
    data_medaka_2 = uq_medaka_2()

    plot_compare(data_rat, data_medaka_1, data_medaka_2)


if __name__ == "__main__":
    compare()

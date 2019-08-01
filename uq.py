from __future__ import absolute_import, division, print_function, unicode_literals

import uncertainpy as un
import chaospy as cp
import numpy as np
import matplotlib.pyplot as plt
import string
import os

from burstiness import is_bursting, is_regular, is_not_spiking, min_spike_amplitude

from uncertainpy.plotting.prettyplot.prettyplot import (
    prettyBar,
    set_latex_font,
    set_style,
    spines_color,
    prettyPlot,
)


# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
figure_format = ".eps"
figure_folder = "figures"


# Uncertainty quantification parameters
polynomial_order = 6
features_to_run = ["is_regular", "is_bursting", "is_not_spiking"]


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
    index = np.arange(1, len(data.uncertain_parameters) + 1) * width

    latex_labels = {
        "g_K": r"$g_\mathrm{K}$",
        "g_Ca": r"$g_\mathrm{Ca}$",
        "g_SK": r"$g_\mathrm{SK}$",
        "g_Na": r"$g_\mathrm{Na}$",
        "g_l": r"$g_\mathrm{l}$",
        "g_BK": r"$g_\mathrm{BK}$",
    }

    xlabels = []
    for label in data.uncertain_parameters:
        xlabels.append(latex_labels[label])

        sensitivity = data[feature].sobol_total_average

        mean = data[feature].mean
        std = np.sqrt(data[feature].variance)
        unit = data[feature].labels

        if std < 1e-14:
            sensitivity = [0] * len(data.uncertain_parameters)

        if len(unit) > 0 and "(" in unit[0]:
            unit = unit[0].split("(")[-1].strip(")")
        else:
            unit = ""

        prettyBar(
            sensitivity,
            xlabels=xlabels,
            nr_colors=len(data.uncertain_parameters),
            index=index,
            ax=ax,
            style=style,
        )

        for tick in ax.get_xticklabels():
            tick.set_rotation(-40)

        ax.set_ylim([0, 1.15])
        # ax.set_title(title, fontsize=titlesize)
        ax.text(
            0.25,
            0.9,
            "Mean = {mean:.2f}".format(mean=mean),
            transform=ax.transAxes,
            fontsize=fontsize,
        )
        ax.text(
            0.25,
            0.8,
            "Std. = {std:.2f}".format(std=std),
            transform=ax.transAxes,
            fontsize=fontsize,
        )

        ax.tick_params(labelsize=fontsize)


def plot(data):
    """
    Plot the total-order Sobol indices.

    Parameters
    ----------
    medaka : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the medaka 1 model.
    """
    style = "seaborn-darkgrid"
    set_style(style)
    set_latex_font()

    # plt.rcParams.update({"axes.titlepad": 12})
    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        squeeze=False,
        sharey="row",
        figsize=(figure_width, figure_width * 0.4),
    )

    axes = axes[0]

    feature_labels = {
        "is_bursting": "IsBursting",
        "is_not_spiking": "IsNotSpiking",
        "is_regular": "IsRegular",
    }

    features = ["is_bursting", "is_regular", "is_not_spiking"]

    yticks = np.arange(0, 1.1, 0.25)

    for i in range(3):
        ax = axes[i]

        plot_sobol_feature(data, features[i], ax)

        ax.text(
            -0.13,
            1.03,
            "\\textbf{{{}}}".format(string.ascii_uppercase[i]),
            transform=ax.transAxes,
            fontsize=titlesize,
        )
        ax.set_title(feature_labels[features[i]])
        ax.set_yticks(yticks)

    axes[0].set_ylabel("Total-order Sobol indices", fontsize=labelsize)

    plt.tight_layout()
    plt.savefig(os.path.join(figure_folder, "sensitivity" + figure_format))

    plt.rcdefaults()


def uq_medaka():
    """
    Perform the uncertainty quantification and sensitivity analysis of the
    medaka model.

    Returns
    -------
    data : uncertainpy.Data
        A data object that contains the results from the uncertainty
        quantification of the rat model.
    """
    parameters = {
        "g_K": 4.18e-4,
        "g_Ca": 6.25e-5,
        "g_SK": 4e-4,
        "g_l": 2e-5,
        "g_BK": 3.13e-4,
        "g_Na": 2.19e-2,
    }
    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))

    # Set full g_BK distribution
    parameters["g_BK"].distribution = cp.Uniform(0, 3.13e-4)

    # Initialize the features
    features = un.SpikingFeatures(
        new_features=[is_bursting, is_regular, is_not_spiking],
        features_to_run=features_to_run,
        logger_level="error",
        strict=False,
        threshold=0.55,
        end_threshold=-0.1,
        normalize=True,
        trim=False,
        min_amplitude=min_spike_amplitude,
    )

    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py", name="medaka", ignore=True)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model, parameters=parameters, features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(
        seed=10,
        filename="medaka",
        polynomial_order=polynomial_order,
        save=True,
        plot=None,
    )

    return data


if __name__ == "__main__":
    data = uq_medaka()

    # To load previous data
    # data = un.Data("data/medaka.h5")
    plot(data)

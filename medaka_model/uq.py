from __future__ import absolute_import, division, print_function, unicode_literals

import uncertainpy as un
import chaospy as cp
import numpy as np
import matplotlib.pyplot as plt
import string

from medaka import scale_conductance
from burstiness import burstiness_factor, bursting, spiking, APs, min_spike_amplitude

from uncertainpy.plotting.prettyplot.prettyplot import prettyBar, set_latex_font, set_style, spines_color, prettyPlot, create_figure


# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
figure_format = ".eps"


# Simulation parameters
simulation_duration = 50000
discard = 10000                 # in ms
simulation_time = simulation_duration + discard
noise_amplitude = 0            # in mV
stimulus_amplitude = 0         # in nA

# Uncertainty quantification parameters
polynomial_order = 5
features_to_run = ["spiking", "bursting", "APs"]



def plot_sobol_feature(data, feature, ax, sobol="total"):
    style = "seaborn-darkgrid"

    width = 0.2
    index = np.arange(1, len(data.uncertain_parameters)+1)*width


    latex_labels = {"g_K": r"$g_\mathrm{K}$",
                    "g_Ca": r"$g_\mathrm{Ca}$",
                    "g_Ca_tabak": r"$g_\mathrm{Ca}$",
                    "g_SK": r"$g_\mathrm{SK}$",
                    "g_Na": r"$g_\mathrm{Na}$",
                    "g_l": r"$g_\mathrm{l}$",
                    "g_BK": r"$g_\mathrm{BK}$"
    }

    xlabels = []
    for label in data.uncertain_parameters:
        xlabels.append(latex_labels[label])

        if feature.lower()  in ["medaka", "medaka_2"]:
            title = "Membrane potential"
        else:
            title = feature.replace("_", " ")
            title = title[0].upper() + title[1:]

        if sobol == "total":
            sensitivity = data[feature].sobol_total_average
        else:
            sensitivity = data[feature].sobol_first_average

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
                  # palette="husl",
                  style=style)

        for tick in ax.get_xticklabels():
            tick.set_rotation(-40)

        ax.set_ylim([0, 1.15])
        ax.set_title(title, fontsize=titlesize)
        ax.text(0.25, 0.9, "Mean = {mean:.2{c}} {unit}".format(mean=mean, c="e" if (abs(mean) < 1e-2 and mean != 0) else "f", unit=unit),
                transform=ax.transAxes, fontsize=fontsize)
        ax.text(0.25, 0.78, "Std. = {std:.2{c}} {unit}".format(std=std, c="e" if (abs(mean) < 1e-2 and mean != 0) else "f", unit=unit),
                transform=ax.transAxes, fontsize=fontsize)

        #ax.set_xticklabels(xlabels, fontsize=labelsize)
        ax.tick_params(labelsize=fontsize)




def plot_compare(tabak, medaka_1, medaka_2, sobol="total"):
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



    features = ["bursting", "spiking", "APs"]
    datas = [tabak, medaka_1, medaka_2]

    yticks = np.arange(0, 1.1, 0.25)

    for i in range(3):
        for j in range(3):
            ax = axes[i][j]

            plot_sobol_feature(datas[i], features[j], ax, sobol=sobol)

            ax.text(-0.17, 1.08, "\\textbf{{{}}}".format(string.ascii_uppercase[i] + str(j + 1)), transform=ax.transAxes, fontsize=titlesize)
            ax.set_title("")
            ax.set_yticks(yticks)


    axes[0][0].set_title("Bursting")
    axes[0][1].set_title("Regular spiking")
    axes[0][2].set_title("AP firing")

    axes[0][0].set_ylabel("RAT", labelpad=12)
    axes[1][0].set_ylabel("MEDAKA 1", labelpad=12)
    axes[2][0].set_ylabel("MEDAKA 2", labelpad=12)

    plt.tight_layout()
    plt.subplots_adjust(left=0.17)
    plt.savefig("sensitivity_" + sobol + figure_format)

    plt.rcdefaults()



def uq_tabak():
    parameters = {"g_K": 9.55e-4,
                  "g_Ca_tabak": 6.37e-4,
                  "g_SK": 6.37e-4,
                  "g_l": 6.37e-5,
                  "g_BK": 2.4e-4}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))


    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, bursting, spiking, APs],
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
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           g_Ca_tabak=6.37e-4,
                           g_Na=0,
                           g_Ca=0)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       plot=None,
                       figure_folder="tabak",
                       filename="tabak",
                       polynomial_order=polynomial_order,
                       save=False)

    return data




def uq_medaka_1():
    parameters = {"g_K": 9.55e-4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": 2.4e-4}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))


    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, bursting, spiking, APs],
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
                           name="medaka",
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
                       figure_folder="medaka_1",
                       filename="medaka_1",
                       polynomial_order=polynomial_order,
                       save=False)

    return data






def uq_medaka_2():
    parameters = {"g_K": 9.55e-4*1.4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4*3,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": 3.2e-4*4*0.67}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their original value
    parameters.set_all_distributions(un.uniform(1))

    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, bursting, spiking, APs],
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

    # # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       figure_folder="medaka_2",
                       filename="medaka_2",
                       plot=None,
                       polynomial_order=polynomial_order,
                       save=False)

    return data




def uq_tabak_full():
    parameters = {"g_K": 9.55e-4,
                  "g_Ca_tabak": 6.37e-4,
                  "g_SK": 6.37e-4,
                  "g_l": 6.37e-5,
                  "g_BK": 2.4e-4}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))

    parameters["g_BK"].distribution = cp.Uniform(0, 3.2e-4)


    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, bursting, spiking, APs],
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
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           g_Ca_tabak=6.37e-4,
                           g_Na=0,
                           g_Ca=0)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       plot=None,
                       figure_folder="tabak",
                       filename="tabak",
                       polynomial_order=polynomial_order,
                       save=False)

    return data


def uq_medaka_1_full():
    parameters = {"g_K": 9.55e-4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": 2.4e-4}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))

    parameters["g_BK"].distribution = cp.Uniform(0, 3.2e-4)


    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, bursting, spiking, APs],
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
                           name="medaka",
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
                       figure_folder="medaka_1",
                       filename="medaka_1",
                       polynomial_order=polynomial_order,
                       save=False)

    return data


def uq_medaka_2_full():
    parameters = {"g_K": 9.55e-4*1.4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4*3,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": 3.2e-4*4*0.67}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))

    parameters["g_BK"].distribution = cp.Uniform(0, 4*3.2e-4)

    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, bursting, spiking, APs],
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

    # # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       figure_folder="medaka_2",
                       filename="medaka_2",
                       plot=None,
                       polynomial_order=polynomial_order,
                       save=False)

    return data




def compare():
    data_tabak = uq_tabak()
    data_medaka_1 = uq_medaka_1()
    data_medaka_2 = uq_medaka_2()

    plot_compare(data_tabak, data_medaka_1, data_medaka_2)
    plot_compare(data_tabak, data_medaka_1, data_medaka_2, sobol="first")


def compare_full():
    data_tabak = uq_tabak_full()
    data_medaka_1 = uq_medaka_1_full()
    data_medaka_2 = uq_medaka_2_full()

    plot_compare(data_tabak, data_medaka_1, data_medaka_2)
    plot_compare(data_tabak, data_medaka_1, data_medaka_2, sobol="first")


if __name__ == "__main__":
    compare()
    # compare_full()

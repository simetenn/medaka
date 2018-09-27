from __future__ import absolute_import, division, print_function, unicode_literals

import uncertainpy as un
import chaospy as cp
import numpy as np
import matplotlib.pyplot as plt
import string

from medaka import scale_conductance
from burstiness import burstiness_factor, burstiness_efel

from uncertainpy.plotting.prettyplot.prettyplot import prettyBar, set_latex_font, set_style, spines_color


# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
figure_format = ".eps"


# Simulation parameters
discard = 1000                 # in ms
simulation_time = 11000        # in ms
noise_amplitude = 0            # in mV


features_to_run = ["burstiness_factor", "spike_rate", "average_AP_overshoot",
                   "average_AHP_depth", "average_duration", "spikiness_factor"]



def spikiness_factor(time, spikes, info):
    """

    Parameters
    ----------
    time : {None, numpy.nan, array_like}
        Time values of the model. If no time values it is None or numpy.nan.
    spikes : uncertainpy.Spikes
        Spikes found in the model result.
    info : dictionary
        Not used, but is required as input

    Returns
    -------
    time : None
    burstiness_factor : int
        If there are any spikes in the voltage trace.


    Calculate the spikiness_factor from a uncertainpy.Spikes object.
    Compatible with uncertainpy.SpikingFeatures.
    """
    flag = 0
    if spikes.nr_spikes > 0:
        flag = 1

    return None, flag




def str_to_latex(text):
    if "_" in text:
        txt = text.split("_")
        return "$" + txt[0] + "_{\mathrm{" + "-".join(txt[1:]) + "}}$"
    else:
        return text


def list_to_latex(texts):
    tmp = []
    for txt in texts:
        tmp.append(str_to_latex(txt))

    return tmp


def plot_sobol(data, data_stim, filename):

    withouth_stim = ["spikiness_factor", "spike_rate"]

    # features = list(data.keys())
    # features.remove("medaka")

    features =  ["spike_rate", "spikiness_factor", "burstiness_factor", "average_duration",
                 "average_AP_overshoot", "average_AHP_depth"]

    nr_plots = len(features)
    grid_size = np.ceil(np.sqrt(nr_plots))
    grid_x_size = int(grid_size)
    grid_y_size = int(np.ceil(nr_plots/float(grid_x_size)))

    set_latex_font()

    style = "seaborn-darkgrid"

    set_style(style)

    set_latex_font()
    plt.rcParams.update({"axes.titlepad": 8})
    fig, axes = plt.subplots(nrows=grid_y_size, ncols=grid_x_size, squeeze=False, sharex='col', sharey='row',
                             figsize=(figure_width, figure_width*0.8))

    labels = data.get_labels("medaka")
    xlabel, ylabel = labels

    # Add a larger subplot to use to set a common xlabel and ylabel
    set_style("seaborn-white")
    ax = fig.add_subplot(111, zorder=-10)
    spines_color(ax, edges={"top": "None", "bottom": "None",
                            "right": "None", "left": "None"})
    ax.tick_params(top='off', bottom='off', left='off', right='off', labelsize=fontsize,
                labelbottom=False, labelleft=False)
    ax.set_ylabel('Total-order Sobol indices', labelpad=30, fontsize=labelsize)


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



    for i in range(0, grid_x_size*grid_y_size):
        nx = i % grid_x_size
        ny = int(np.floor(i/float(grid_x_size)))

        ax = axes[ny][nx]

        if i < nr_plots:
            if features[i] == "medaka":
                title = "Membrane potential"
            else:
                title = features[i].replace("_", " ")
                title = title[0].upper() + title[1:]

            if features[i] in withouth_stim:
                sensitivity = data[features[i]].sobol_total_average
                mean = data[features[i]].mean
                std = np.sqrt(data[features[i]].variance)
                unit = data[features[i]].labels
            else:
                sensitivity = data_stim[features[i]].sobol_total_average
                mean = data_stim[features[i]].mean
                std = np.sqrt(data_stim[features[i]].variance)
                unit = data_stim[features[i]].labels


            if len(unit) > 0 and "(" in unit[0]:
                unit = unit[0].split("(")[-1].strip(")")
            else:
                unit = ""

            if features[i] == "spike_rate":
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
            # ax.set_title("({}) ".format(string.ascii_uppercase[i]) + title, fontsize=titlesize)
            ax.set_title(title, fontsize=titlesize)
            ax.text(-0.08, 1.08, "\\textbf{{{}}}".format(string.ascii_uppercase[i]), transform=ax.transAxes, fontsize=titlesize)

            ax.text(0.2, 0.9, "Mean = {mean:.2{c}} {unit}".format(mean=mean, c="e" if abs(mean) < 1e-2 else "f", unit=unit),
                    transform=ax.transAxes, fontsize=labelsize)
            ax.text(0.2, 0.8, "Std. = {std:.2{c}} {unit}".format(std=std, c="e" if abs(mean) < 1e-2 else "f", unit=unit),
                    transform=ax.transAxes, fontsize=labelsize)

            #ax.set_xticklabels(xlabels, fontsize=labelsize)
            ax.tick_params(labelsize=fontsize)
        else:
            ax.axis("off")

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.13, left=0.1)
    plt.savefig(filename + figure_format)

    plt.rcdefaults()



def burstiness_robustness():
    parameters = {"g_K": scale_conductance(3),
                  "g_Ca": scale_conductance(2),
                  "g_SK": scale_conductance(2),
                  "g_l": scale_conductance(0.2)}

    # Create the parameters
    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    # within a +/- 50% interval around their fixed value
    parameters.set_all_distributions(un.uniform(1))


    # Initialize the features
    features = un.SpikingFeatures(new_features=burstiness_factor, strict=False,
                                  logger_level="error", features_to_run="burstiness_factor")



    # g_bk => 0
    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py", name="medaka",
                           discard=discard, noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           gbk_bk=scale_conductance(0))

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data_0 = UQ.quantify(seed=10, plot=None)


    # g_bk => 0.5
    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py", name="medaka",
                           discard=discard, noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           gbk_bk=scale_conductance(0.5))

    # Switch to the model with new default values
    UQ.model = model

    # We set the seed to easier be able to reproduce the result
    data_05 = UQ.quantify(seed=10, plot=None)


    # g_bk => 0
    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py", name="medaka",
                           discard=discard, noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           gbk_bk=scale_conductance(1))

    # Switch to the model with new default values
    UQ.model = model

    # We set the seed to easier be able to reproduce the result
    data_1 = UQ.quantify(seed=10, plot=None)


    # plot the results
    mean = [data_0["burstiness_factor"].mean,
            data_05["burstiness_factor"].mean,
            data_1["burstiness_factor"].mean]

    error = [np.sqrt(data_0["burstiness_factor"].variance),
             np.sqrt(data_05["burstiness_factor"].variance),
             np.sqrt(data_1["burstiness_factor"].variance)]

    set_latex_font()

    xlabels = [r"$g_{BK} = 0$ (mS/cm$^2$)" + "\n" + r"$\rightarrow 0$ nS",
               r"$g_{BK} = 0.16$ (mS/cm$^2$)" + "\n" + r"$\rightarrow 0.5$ nS",
               r"$g_{BK} = 0.32$ (mS/cm$^2$)" + "\n" + r"$\rightarrow 1$ nS"]

    ax = prettyBar(mean, error, xlabels=xlabels, style="seaborn-darkgrid", color=0, width=0.5, ylabel="Burstiness")

    plt.savefig("burstiness_robustness" + figure_format)













def burstiness_tabak_vs_medaka():

    # Tabak
    parameters = {"g_K": 9.55e-4,
                  "g_Ca_tabak": 6.37e-4,
                  "g_SK": 6.37e-4,
                  "g_l": 6.37e-5,
                  "g_BK": scale_conductance(0.67)}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))

    # Initialize the features
    features = un.SpikingFeatures(new_features=burstiness_factor,
                                  features_to_run="burstiness_factor",
                                  logger_level="error",
                                  strict=False,
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False)


    # Initialize the Tabak model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           g_Na=0,
                           g_Ca=0,
                           stimulus_amplitude=0.0015)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)
    data_tabak = UQ.quantify(seed=10,
                             plot="all",
                             figure_folder="tabak",
                             filename="tabak",
                             polynomial_order=5)


    # Medaka
    parameters = {"g_K": 9.55e-4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": scale_conductance(0.67)}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))


    # Initialize the Medaka model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           stimulus_amplitude=0.0015)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data_medaka= UQ.quantify(seed=10,
                             figure_folder="medaka_1_stim",
                             filename="medaka_1_stim",
                             plot="all",
                             polynomial_order=5)




    # Plotting
    style = "seaborn-darkgrid"
    set_style(style)
    set_latex_font()

    plt.rcParams.update({"axes.titlepad": 8})
    fig, axes = plt.subplots(nrows=1, ncols=2, squeeze=False, sharex="col", sharey="row",
                             figsize=(figure_width, figure_width*0.4))

    ax1 = axes[0][0]
    ax2 = axes[0][1]

    ax1.set_ylabel("Total-order Sobol indices", labelpad=30, fontsize=labelsize)


    width = 0.2

    latex_labels = {"g_K": r"$g_\mathrm{K}$",
                    "g_Ca": r"$g_\mathrm{Ca}$",
                    "g_Ca_tabak": r"$g_\mathrm{Ca}$",
                    "g_SK": r"$g_\mathrm{SK}$",
                    "g_Na": r"$g_\mathrm{Na}$",
                    "g_l": r"$g_\mathrm{l}$",
                    "g_BK": r"$g_\mathrm{BK}$"
    }




    # Tabak
    xlabels = []
    for label in data_tabak.uncertain_parameters:
        xlabels.append(latex_labels[label])

    index = np.arange(1, len(data_tabak.uncertain_parameters)+1)*width


    sensitivity = data_tabak["burstiness_factor"].sobol_total_average
    mean = data_tabak["burstiness_factor"].mean
    std = np.sqrt(data_tabak["burstiness_factor"].variance)

    prettyBar(sensitivity,
              xlabels=xlabels,
              title="Tabak",
              nr_colors=len(data_tabak.uncertain_parameters),
              index=index,
              ax=ax1,
              style=style)

    for tick in ax1.get_xticklabels():
        tick.set_rotation(-40)

    ax1.set_ylim([0, 1])
    # ax1.set_title(title, fontsize=titlesize)
    ax1.text(-0.08, 1.08, R"\textbf{A}", transform=ax1.transAxes, fontsize=titlesize)

    ax1.text(0.38, 0.9, "Mean = {mean:.2{c}}".format(mean=mean, c="e" if abs(mean) < 1e-2 else "f"),
            transform=ax1.transAxes, fontsize=labelsize)
    ax1.text(0.38, 0.8, "Std. = {std:.2{c}}".format(std=std, c="e" if abs(mean) < 1e-2 else "f"),
            transform=ax1.transAxes, fontsize=labelsize)

    ax1.tick_params(labelsize=fontsize)



    # Medaka
    xlabels = []
    for label in data_medaka.uncertain_parameters:
        xlabels.append(latex_labels[label])

    index = np.arange(1, len(data_medaka.uncertain_parameters)+1)*width


    sensitivity = data_medaka["burstiness_factor"].sobol_total_average
    mean = data_medaka["burstiness_factor"].mean
    std = np.sqrt(data_medaka["burstiness_factor"].variance)

    prettyBar(sensitivity,
              xlabels=xlabels,
              title="Medaka",
              nr_colors=len(data_medaka.uncertain_parameters),
              index=index,
              ax=ax2,
              style=style)

    for tick in ax2.get_xticklabels():
        tick.set_rotation(-40)

    ax2.set_ylim([0, 1.15])
    ax2.text(-0.08, 1.08, R"\textbf{B}", transform=ax2.transAxes, fontsize=titlesize)

    ax2.text(0.38, 0.9, "Mean = {mean:.2{c}}".format(mean=mean, c="e" if abs(mean) < 1e-2 else "f"),
            transform=ax2.transAxes, fontsize=labelsize)
    ax2.text(0.38, 0.8, "Std. = {std:.2{c}}".format(std=std, c="e" if abs(mean) < 1e-2 else "f"),
            transform=ax2.transAxes, fontsize=labelsize)

    ax2.tick_params(labelsize=fontsize)


    plt.tight_layout()
    plt.savefig("tabak_vs_medaka" + figure_format)

    plt.rcdefaults()




def uncertain_tabak():
    parameters = {"g_K": 9.55e-4,
                  "g_Ca_tabak": 6.37e-4,
                  "g_SK": 6.37e-4,
                  "g_l": 6.37e-5,
                  "g_BK": scale_conductance(0.67)}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))

    # parameters["g_Na"].distribution = cp.Uniform(0.07,
    #                                              0.07+0.07/2)

    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, spikiness_factor],
                                  features_to_run=features_to_run,
                                  logger_level="error",
                                  strict=False,
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False)


    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           g_Na=0,
                           g_Ca=0,
                           stimulus_amplitude=0.0015)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       plot="all",
                       figure_folder="tabak",
                       filename="tabak",
                       polynomial_order=5)


    plot_sobol(data, data, "sensitivity_tabak")








def uncertain_medaka():
    parameters = {"g_K": 9.55e-4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": scale_conductance(0.67)}

    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))

    # parameters["g_Na"].distribution = cp.Uniform(0.07,
    #                                              0.07+0.07/2)

    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, spikiness_factor],
                                  features_to_run=features_to_run,
                                  logger_level="error",
                                  strict=False,
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False)

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
                       plot="all",
                       figure_folder="medaka_1",
                       filename="medaka_1",
                       polynomial_order=5)



    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           stimulus_amplitude=0.0015)

    # Perform the uncertainty quantification
    UQ.model = model

    # We set the seed to easier be able to reproduce the result
    data_stim = UQ.quantify(seed=10,
                            figure_folder="medaka_1_stim",
                            filename="medaka_1_stim",
                            plot="all",
                            polynomial_order=5)

    plot_sobol(data, data_stim, "sensitivity_medaka_1")






def uncertain_medaka_2():
    parameters = {"g_K": 9.55e-4*1.4,
                  "g_Ca": 2e-4,
                  "g_SK": 6.37e-4*3,
                  "g_Na": 0.07,
                  "g_l": 6.37e-5,
                  "g_BK": scale_conductance(0.67)*3}


    parameters = un.Parameters(parameters)

    # Set all parameters to have a uniform distribution
    parameters.set_all_distributions(un.uniform(1))

    # Initialize the features
    features = un.SpikingFeatures(new_features=[burstiness_factor, spikiness_factor],
                                  features_to_run=features_to_run,
                                  strict=False,
                                  logger_level="error",
                                  threshold=0.55,
                                  end_threshold=-0.1,
                                  normalize=True,
                                  trim=False)

    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           taun=5,
                           vf=-15)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10,
                       figure_folder="medaka_2",
                       filename="medaka_2",
                       plot="all",
                       polynomial_order=5)

    # Initialize the model and defining default options
    model = un.NeuronModel(file="medaka.py",
                           name="medaka",
                           discard=discard,
                           noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           ignore=True,
                           taun=5,
                           vf=-15,
                           stimulus_amplitude=0.0015)

    # Perform the uncertainty quantification
    UQ.model = model

    # We set the seed to easier be able to reproduce the result
    data_stim = UQ.quantify(seed=10,
                            figure_folder="medaka_2_stim",
                            filename="medaka_2_stim",
                            plot="all",
                            polynomial_order=5)

    plot_sobol(data, data_stim, "sensitivity_medaka_2")




if __name__ == "__main__":
    # burstiness_robustness()
    burstiness_tabak_vs_medaka()
    # uncertain_tabak()
    # uncertain_medaka()
    uncertain_medaka_2()
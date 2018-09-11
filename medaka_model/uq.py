from __future__ import absolute_import, division, print_function, unicode_literals

import uncertainpy as un
import chaospy as cp
import numpy as np
import matplotlib.pyplot as plt

from medaka import scale_conductance
from burstiness import burstiness_factor, burstiness_efel

from uncertainpy.plotting.prettyplot.prettyplot import prettyBar, set_latex_font


# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
figure_format = ".png"


# Simulation parameters
discard = 1000                 # in ms
simulation_time = 11000        # in ms
noise_amplitude = 0            # in mV


def burstiness_robustness():
    parameters = {"gkdrbar_kdrt": scale_conductance(3.2),
                  "ghvat_ihvat": scale_conductance(2),
                  "gskbar_sk": scale_conductance(2),
                  "g_pas": scale_conductance(0.2)}

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
    model = un.NeuronModel(path="", file="medaka.py", name="medaka",
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
    model = un.NeuronModel(path="", file="medaka.py", name="medaka",
                           discard=discard, noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time,
                           gbk_bk=scale_conductance(0.5))

    # Switch to the model with new default values
    UQ.model = model

    # We set the seed to easier be able to reproduce the result
    data_05 = UQ.quantify(seed=10, plot=None)


    # g_bk => 0
    # Initialize the model and defining default options
    model = un.NeuronModel(path="", file="medaka.py", name="medaka",
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



def uncertain_gbk():
    parameters = {"gkdrbar_kdrt": scale_conductance(3.2),
                  "ghvat_ihvat": scale_conductance(2),
                  "gskbar_sk": scale_conductance(2),
                  "g_pas": scale_conductance(0.2),
                  "gbk_bk": scale_conductance(1)}

    parameters = un.Parameters(parameters)

    # # Set all parameters to have a uniform distribution
    # # within a 20% interval around their fixed value
    parameters.set_all_distributions(un.uniform(1))

    # parameters["gbk_bk"].distribution = cp.Uniform(scale_conductance(0),
    #                                                scale_conductance(1))


    # parameters = {"gkdrbar_kdrt": cp.Uniform(scale_conductance(2.8),
    #                                          scale_conductance(3.8)),
    #               "ghvat_ihvat": cp.Uniform(scale_conductance(1.5),
    #                                          scale_conductance(2.5)),
    #               "gskbar_sk": cp.Uniform(scale_conductance(1.4),
    #                                       scale_conductance(2.5)),
    #               "g_pas": cp.Uniform(scale_conductance(0),
    #                                   scale_conductance(1)),
    #               "gbk_bk": cp.Uniform(scale_conductance(0),
    #                                    scale_conductance(1))}

    # Initialize the features
    features = un.SpikingFeatures(new_features=burstiness_factor, strict=False,
                                  logger_level="error", features_to_run="all")

    # efel_features_to_run = ["burst_mean_freq", "AP_width", "spike_width2", "burstiness_efel"]

    # features = un.EfelFeatures(features_to_run=efel_features_to_run,
    #                            new_features=burstiness_efel,
    #                            strict=False, logger_level="error")

    # Initialize the model and defining default options
    model = un.NeuronModel(path="", file="medaka.py", name="medaka",
                           discard=discard, noise_amplitude=noise_amplitude,
                           simulation_time=simulation_time)

    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)

    # We set the seed to easier be able to reproduce the result
    data = UQ.quantify(seed=10, plot="all")


if __name__ == "__main__":
    # burstiness_robustness()
    uncertain_gbk()
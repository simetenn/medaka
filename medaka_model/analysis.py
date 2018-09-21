from __future__ import absolute_import, division, print_function, unicode_literals

import uncertainpy as un
import chaospy as cp
import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib.patches import Patch

from medaka import scale_conductance, medaka
from burstiness import burstiness, duration, burst_threshold


data_folder = "../data"

# Simulation and analysis parameters
discard = 1000                 # in ms
simulation_time = 50000        # in ms
noise_amplitude = 0.004        # in nA
robustness_reruns = 512        # From original article
simulation_time_plot = 5000    # in ms


# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
figure_format = ".eps"
label_x = -0.08
label_y = 1.08
axis_grey = (0.6, 0.6, 0.6)


# Set the random seed to increase reproducability
np.random.seed(10)


# Set default options for plotting
params = {
    "text.usetex": True,
    "text.latex.preamble": "\\usepackage{lmodern}, \\usepackage{amsmath}, \\usepackage{amssymb}",
    "font.family": "lmodern",
    "xtick.color": axis_grey,
    "ytick.color": axis_grey,
    "axes.edgecolor": axis_grey,
    "xtick.bottom": True,
    "ytick.left": True,
    "axes.spines.bottom": True,
    "axes.spines.left": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1,
    "axes.labelsize": labelsize
}


def calculate_frequency_bf(**parameters):
    """
    Calculate the frequency of a binned burstiness factor from a series of
    model evaluations.

    Parameters
    ----------
    **parameters
        Any number of optional parameters passed on to the model.

    Returns
    -------
    bins : numpy.array
        The bins for the event durations.
    frequency : numpy.array
        Frequency of occurrences in each binned event duration
    burstiness_factor : float
        The fraction of events that is above `burst_threshold`, for all
        reruns.

    Notes
    -----
    """
    nr_bins = 40
    hist_range = (0, 250)
    time, voltage = medaka(noise_amplitude=noise_amplitude,
                           discard=discard,
                           simulation_time=simulation_time,
                           **parameters)

    event_durations = duration(time, voltage)

    binned_durations, bins = np.histogram(event_durations, bins=nr_bins, range=hist_range)
    frequency = binned_durations/binned_durations.sum()

    burstiness_factor = burstiness(event_durations)

    return bins, frequency, burstiness_factor



def change_g_BK(**parameters):
    """
    Change g_BK values and calculate the burstiness factor for each.

    Parameters
    ----------
    **parameters
        Parameters for the model, can not be g_BK.

    Returns
    -------
    g_BKs : numpy.array
        The g_BK values in S/cm^2.
    burstiness_factors : list
        The burstiness factor for each g_BK.
    """
    original_g_BKs = np.arange(0, 1.1, 0.1)

    g_BKs = scale_conductance(original_g_BKs)

    burstiness_factors = []

    for g_BK in g_BKs:
        bins, frequency, burstiness_factor = calculate_frequency_bf(g_BK=g_BK, **parameters)
        burstiness_factors.append(burstiness_factor)

    return g_BKs, burstiness_factors



def change_tau_BK(**parameters):
    """
    Change tau_BK values and calculate the burstiness factor for each.

    Parameters
    ----------
    **parameters
        Parameters for the model, can not be g_BK or tau_BK.

    Returns
    -------
    tau_BKs : numpy.array
        The g_BK values in S/cm^2.
    burstiness_factors : list
        The burstiness factor for each tau_BK.

    Notes
    -----
    Uses original g_BK => 1.
    """
    tau_BKs = np.arange(2, 10.1, 1)

    burstiness_factors = []

    for tau_BK in tau_BKs:
        bins, frequency, burstiness_factor = calculate_frequency_bf(tau_BK=tau_BK, **parameters)
        burstiness_factors.append(burstiness_factor)

    return tau_BKs, burstiness_factors



def robustness(g_BK=0):
    """
    Calculate the number of occurrences for binned burstiness factor of several
    model runs with varying conductances (except g_BK)

    Parameters
    ----------
    g_BKs : float
        The value of the g_BK conductance in S/cm^2.


    Returns
    -------
        Returns
    -------
    bins : numpy.array
        The bins for the burstiness.
    binned_burstiness_factors : numpy.array
        The number of model evaluations with burstiness factor corresponding to
        each bin.
    """
    bins = 10
    hist_range = (0, 1)

    # Original values (scaled to the new model)
    g_l_scaled = scale_conductance(0.2)
    g_K_scaled = scale_conductance(3)
    g_Ca_scaled = scale_conductance(2)
    g_SK_scaled = scale_conductance(2)

    # Draw conductances from uniform distributions +/- 50% of their original values
    g_K = np.random.uniform(g_K_scaled*0.5, g_K_scaled*1.5, robustness_reruns)
    g_Ca = np.random.uniform(g_Ca_scaled*0.5, g_Ca_scaled*1.5, robustness_reruns)
    g_SK = np.random.uniform(g_SK_scaled*0.5, g_SK_scaled*1.5, robustness_reruns)
    g_l = np.random.uniform(g_l_scaled*0.5, g_l_scaled*1.5, robustness_reruns)

    burstiness_factors = []

    # Run the model for each of the selected conductances
    # and calculate the burstiness factor of each evaluation
    for i in range(robustness_reruns):
        time, voltage = medaka(noise_amplitude=noise_amplitude,
                               discard=discard,
                               simulation_time=simulation_time,
                               g_BK=g_BK,
                               g_K=g_K[i],
                               g_Ca=g_Ca[i],
                               g_SK=g_SK[i],
                               g_l=g_l[i])

        event_durations = duration(time, voltage)

        burstiness_factor = burstiness(event_durations)

        if burstiness_factor is not None:
            burstiness_factors.append(burstiness_factor)

    binned_burstiness_factors, bins = np.histogram(burstiness_factors, bins=bins, range=hist_range)

    return bins, binned_burstiness_factors




def figure_1():
    """
    Recreate figure 1 in Tabak et. al. 2011. Figure is saved as figure_1.png

    http://www.jneurosci.org/content/31/46/16855/tab-article-info
    """
    # g_bk => 0
    g_BK = scale_conductance(0)

    # Medaka
    time_0, V_0 = medaka(noise_amplitude=noise_amplitude,
                         discard=discard,
                         g_BK=g_BK,
                         simulation_time=simulation_time)

    bins_0, frequency_0, burstiness_factor_0 = calculate_frequency_bf(g_BK=g_BK)

    # Tabak
    time_0_tabak, V_0_tabak = medaka(noise_amplitude=noise_amplitude,
                         discard=discard,
                         g_BK=g_BK,
                         g_Ca_tabak=6.37e-4,
                         g_Na=0,
                         g_Ca=0,
                         simulation_time=simulation_time)

    bins_0_tabak, frequency_0_tabak, burstiness_factor_0_tabak = calculate_frequency_bf(g_BK=g_BK,
                                                                                        g_Ca_tabak=6.37e-4,
                                                                                        g_Na=0,
                                                                                        g_Ca=0)



    # g_bk => 0.5
    g_BK = scale_conductance(0.5)

    # Medaka
    time_05, V_05 = medaka(noise_amplitude=noise_amplitude,
                           discard=discard,
                           g_BK=g_BK,
                           simulation_time=simulation_time)

    bins_05, frequency_05, burstiness_factor_05 = calculate_frequency_bf(g_BK=g_BK)

    # Tabak
    time_05_tabak, V_05_tabak = medaka(noise_amplitude=noise_amplitude,
                         discard=discard,
                         g_BK=g_BK,
                         g_Ca_tabak=6.37e-4,
                         g_Na=0,
                         g_Ca=0,
                         simulation_time=simulation_time)

    bins_05_tabak, frequency_05_tabak, burstiness_factor_05_tabak = calculate_frequency_bf(g_BK=g_BK,
                                                                                           g_Ca_tabak=6.37e-4,
                                                                                           g_Na=0,
                                                                                           g_Ca=0)


    # g_bk => 0.65
    g_BK = scale_conductance(0.7)
    # Medak
    time_065, V_065 = medaka(noise_amplitude=noise_amplitude,
                           discard=discard,
                           g_BK=g_BK,
                           simulation_time=simulation_time)

    bins_065, frequency_065, burstiness_factor_065 = calculate_frequency_bf(g_BK=g_BK)

    # Tabak
    time_065_tabak, V_065_tabak = medaka(noise_amplitude=noise_amplitude,
                         discard=discard,
                         g_BK=g_BK,
                         g_Ca_tabak=6.37e-4,
                         g_Na=0,
                         g_Ca=0,
                         simulation_time=simulation_time)

    bins_065_tabak, frequency_065_tabak, burstiness_factor_065_tabak = calculate_frequency_bf(g_BK=g_BK,
                                                                                              g_Ca_tabak=6.37e-4,
                                                                                              g_Na=0,
                                                                                              g_Ca=0)


    # g_bk => 1
    g_BK = scale_conductance(1)
    # Medaka
    time_1, V_1 = medaka(noise_amplitude=noise_amplitude,
                         discard=discard,
                         g_BK=g_BK,
                         simulation_time=simulation_time)

    bins_1, frequency_1, burstiness_factor_1 = calculate_frequency_bf(g_BK=g_BK)

    # Tabak
    time_1_tabak, V_1_tabak = medaka(noise_amplitude=noise_amplitude,
                         discard=discard,
                         g_BK=g_BK,
                         g_Ca_tabak=6.37e-4,
                         g_Na=0,
                         simulation_time=simulation_time)

    bins_1_tabak, frequency_1_tabak, burstiness_factor_1_tabak = calculate_frequency_bf(g_BK=g_BK,
                                                                                        g_Ca_tabak=6.37e-4,
                                                                                        g_Na=0,
                                                                                        g_Ca=0)

    # Calculate results for figure 1D
    # Medaka
    scaled_g_BKs, burstiness_factors_g_BK = change_g_BK()
    # Tabak
    scaled_g_BKs_tabak, burstiness_factors_g_BK_tabak = change_g_BK(g_Ca_tabak=6.37e-4,
                                                                    g_Na=0,
                                                                    g_Ca=0)



    # Calculate results for figure 1E
    # Medaka
    scaled_tau_BK, burstiness_factors_tau_BK = change_tau_BK(g_BK=scale_conductance(0))

    # Tabak
    scaled_tau_BK_tabak, burstiness_factors_tau_BK_tabak = change_tau_BK(g_Ca_tabak=6.37e-4,
                                                                         g_Na=0,
                                                                         g_Ca=0,
                                                                         g_BK=scale_conductance(1))


    # Plotting
    plt.rcParams.update(params)

    fig = plt.figure(figsize=(figure_width, figure_width))
    gs = plt.GridSpec(4, 6)
    ax1 = plt.subplot(gs[0, :-2])
    ax2 = plt.subplot(gs[1, :-2])
    ax3 = plt.subplot(gs[2, :-2])

    ax4 = plt.subplot(gs[0, -2:])
    ax5 = plt.subplot(gs[1, -2:])
    ax6 = plt.subplot(gs[2, -2:])

    ax7 = plt.subplot(gs[3, :3])
    ax8 = plt.subplot(gs[3, 3:])


    voltage_axes = [ax1, ax2, ax3]
    burst_axes = [ax4, ax5, ax6]



    ax1.plot(time_0_tabak, V_0_tabak, color="tab:gray")
    ax1.plot(time_0, V_0, color="tab:blue")
    title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(0)*1000) + r"$ (mS/cm$^2$)"
    ax1.set_title(title)
    ax1.get_xaxis().set_visible(False)
    ax1.text(label_x, label_y, r"\textbf{A}", transform=ax1.transAxes, fontsize=titlesize)

    ax2.plot(time_05_tabak, V_05_tabak, color="tab:gray")
    ax2.plot(time_05, V_05, color="tab:blue")
    title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(0.5)*1000) + r"$ (mS/cm$^2$) "
    ax2.set_title(title)
    ax2.get_xaxis().set_visible(False)
    ax2.text(label_x, label_y, r"\textbf{B}", transform=ax2.transAxes, fontsize=titlesize)


    ax3.plot(time_1_tabak, V_1_tabak, color="tab:gray")
    ax3.plot(time_1, V_1, color="tab:blue")
    title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(1)*1000) + r"$ (mS/cm$^2$)"
    ax3.set_title(title)
    ax3.set_xlabel("Time (ms)", fontsize=labelsize)
    ax3.text(label_x, label_y, r"\textbf{C}", transform=ax3.transAxes, fontsize=titlesize)


    yticks = [-60, -40, -20, 0, 20]

    for ax in voltage_axes:
        ax.set_ylabel("V (mV)")
        ax.set_ylim([-75, 40])
        ax.set_xlim([discard, simulation_time_plot])
        ax.set_yticks(yticks)
        ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")



    ax4.bar(bins_0_tabak[:-1],
            frequency_0_tabak,
            width=(bins_0_tabak[1] - bins_0_tabak[0]),
            align="edge",
            color="tab:grey")
    ax4.text(100, 0.85, "BF = {:.2f}".format(burstiness_factor_0_tabak), color="tab:grey")

    ax4.bar(bins_0[:-1],
            frequency_0,
            width=(bins_0[1] - bins_0[0]),
            align="edge",
            color="tab:blue")
    ax4.text(100, 0.7, "BF = {:.2f}".format(burstiness_factor_0), color="tab:blue")


    ax5.bar(bins_05[:-1],
            frequency_05,
            width=(bins_05[1] - bins_05[0]),
            align="edge",
            color="tab:blue")
    ax5.text(100, 0.7, "BF = {:.2f}".format(burstiness_factor_05), color="tab:blue")

    ax5.bar(bins_05_tabak[:-1],
            frequency_05_tabak,
            width=(bins_05_tabak[1] - bins_05_tabak[0]),
            align="edge",
            color="tab:grey")
    ax5.text(100, 0.85, "BF = {:.2f}".format(burstiness_factor_05_tabak), color="tab:grey")


    ax6.bar(bins_1_tabak[:-1],
            frequency_1_tabak,
            width=(bins_1_tabak[1] - bins_1_tabak[0]),
            align="edge",
            color="tab:grey")
    ax6.text(100, 0.85, "BF = {:.2f}".format(burstiness_factor_1_tabak), color="tab:grey")

    ax6.bar(bins_1[:-1],
            frequency_1,
            width=(bins_1[1] - bins_1[0]),
            align="edge")
    ax6.text(100, 0.7, "BF = {:.2f}".format(burstiness_factor_1), color="tab:blue")


    yticks = [0, 0.2,  0.4,  0.6,  0.8, 1]
    # xticks = [0, 0.05,  0.1,  0.15,  0.2, 0.25]
    xticks = [0, 50,  100,  150,  200, 250]


    for ax in burst_axes:
        ax.axvline(burst_threshold, color=axis_grey)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, .3])
        ax.set_yticks(yticks)
        ax.set_xticks(xticks)
        ax.set_ylabel("Frequency")
        ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")

    ax6.set_xlabel("Event duration (ms)")


    ax7.plot(scaled_g_BKs_tabak*1000, burstiness_factors_g_BK_tabak, marker=".", color="tab:grey")
    ax7.plot(scaled_g_BKs*1000, burstiness_factors_g_BK, marker=".", color="tab:blue")
    ax7.set_xlabel(r"$g_{BK}$ (mS/cm$^2$)")
    ax7.set_ylabel("Burstiness")
    ax7.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")
    ax7.set_yticks(yticks)
    ax7.set_ylim([-0.05, 1.05])
    # ax7.set_xlim([0, scale_conductance(1)*1000])
    # ax7.set_xlabel(r"$g_{BK}$ (mS/cm$^2$)")


    ax8.plot(scaled_tau_BK_tabak, burstiness_factors_tau_BK_tabak, marker=".", color="tab:grey")
    ax8.plot(scaled_tau_BK, burstiness_factors_tau_BK, marker=".", color="tab:blue")
    ax8.set_xlabel(r"$\tau_{BK}$ (ms)")
    ax8.set_ylabel("Burstiness")
    ax8.set_yticks(yticks)
    ax8.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")
    ax8.set_ylim([-0.05, 1.05])
    # ax8.set_xlim([2, 20])

    ax7.text(label_x, label_y, r"\textbf{D}", transform=ax7.transAxes, fontsize=titlesize)
    ax8.text(label_x, label_y, r"\textbf{E}", transform=ax8.transAxes, fontsize=titlesize)


    rat_legend = Patch(facecolor="tab:grey", label="Rat")
    fish_legend = Patch(facecolor="tab:blue", label="Medaka 1")
    plt.legend(handles=[rat_legend, fish_legend],
               bbox_to_anchor=(0.65, 1),
               bbox_transform=plt.gcf().transFigure,
               ncol=2)


    plt.tight_layout()

    plt.subplots_adjust(top=0.91)

    plt.savefig("figure_1" + figure_format)




# def figure_1():
#     """
#     Recreate figure 1 in Tabak et. al. 2011. Figure is saved as figure_1.png

#     http://www.jneurosci.org/content/31/46/16855/tab-article-info
#     """
#     simulation_time_plot = 5000

#     ## g_bk => 0
#     g_BK = scale_conductance(0)

#     # Medaka
#     time_0, V_0 = medaka(noise_amplitude=noise_amplitude,
#                          discard=discard,
#                          g_BK=g_BK,
#                          simulation_time=simulation_time)

#     event_durations_0 = duration(time_0, V_0)

#     bins_0, frequency_0, burstiness_factor_0 = calculate_frequency_bf(g_BK=g_BK)

#     # Tabak

#     ## g_bk => 0.5
#     g_BK = scale_conductance(0.5)
#     time_05, V_05 = medaka(noise_amplitude=noise_amplitude,
#                            discard=discard,
#                            g_BK=g_BK,
#                            simulation_time=simulation_time)

#     event_durations_05 = duration(time_05, V_05)

#     bins_05, frequency_05, burstiness_factor_05 = calculate_frequency_bf(g_BK=g_BK)



#     # g_bk => 1
#     g_BK = scale_conductance(1)
#     time_1, V_1 = medaka(noise_amplitude=noise_amplitude,
#                          discard=discard,
#                          g_BK=g_BK,
#                          simulation_time=simulation_time)

#     event_durations_1 = duration(time_1, V_1)

#     bins_1, frequency_1, burstiness_factor_1 = calculate_frequency_bf(g_BK=g_BK)



#     # Calculate results for figure 1D
#     scaled_g_BKs, burstiness_factors_g_BK = change_g_BK()


#     # Calculate results for figure 1E
#     scaled_tau_BK, burstiness_factors_tau_BK = change_tau_BK()


#     # Rescale from ms to s
#     time_0 /= 1000
#     time_05 /= 1000
#     time_1 /= 1000
#     bins_0 /= 1000
#     bins_05 /= 1000
#     bins_1 /= 1000
#     burst_threshold_scaled = burst_threshold/1000
#     simulation_time_plot_scaled = simulation_time_plot/1000
#     discard_scaled = discard/1000


#     # Plotting
#     plt.rcParams.update(params)

#     fig = plt.figure(figsize=(figure_width, figure_width))
#     gs = plt.GridSpec(4, 6)
#     ax1 = plt.subplot(gs[0, :-2])
#     ax2 = plt.subplot(gs[1, :-2])
#     ax3 = plt.subplot(gs[2, :-2])

#     ax4 = plt.subplot(gs[0, -2:])
#     ax5 = plt.subplot(gs[1, -2:])
#     ax6 = plt.subplot(gs[2, -2:])

#     ax7 = plt.subplot(gs[3, :3])
#     ax8 = plt.subplot(gs[3, 3:])


#     voltage_axes = [ax1, ax2, ax3]
#     burst_axes = [ax4, ax5, ax6]

#     ax1.plot(time_0, V_0)
#     title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(0)*1000) + r"$ (mS/cm$^2$) $\rightarrow 0$ nS"
#     ax1.set_title(title)
#     ax1.get_xaxis().set_visible(False)
#     ax1.text(label_x, label_y, r"\textbf{A}", transform=ax1.transAxes, fontsize=titlesize)

#     ax2.plot(time_05, V_05)
#     title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(0.5)*1000) + r"$ (mS/cm$^2$) $\rightarrow 0.5$ nS"
#     ax2.set_title(title)
#     ax2.get_xaxis().set_visible(False)
#     ax2.text(label_x, label_y, r"\textbf{B}", transform=ax2.transAxes, fontsize=titlesize)

#     ax3.plot(time_1, V_1)
#     title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(1)*1000) + r"$ (mS/cm$^2$) $\rightarrow 1$ nS"
#     ax3.set_title(title)
#     ax3.set_xlabel("Time (s)", fontsize=labelsize)
#     ax3.text(label_x, label_y, r"\textbf{C}", transform=ax3.transAxes, fontsize=titlesize)


#     yticks = [-60, -40, -20, 0]

#     for ax in voltage_axes:
#         ax.set_ylabel("V (mV)")
#         ax.set_ylim([-70, 10])
#         ax.set_xlim([discard_scaled, simulation_time_plot_scaled])
#         ax.set_yticks(yticks)
#         ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")



#     ax4.bar(bins_0[:-1], frequency_0, width=(bins_0[1] - bins_0[0]), align="edge")
#     ax4.text(0.1, 0.8, "BF = {}".format(burstiness_factor_0))

#     ax5.bar(bins_05[:-1], frequency_05, width=(bins_05[1] - bins_05[0]), align="edge")
#     ax5.text(0.1, 0.8, "BF = {:.2f}".format(burstiness_factor_05))

#     ax6.bar(bins_1[:-1], frequency_1, width=(bins_1[1] - bins_1[0]), align="edge")
#     ax6.text(0.1, 0.8, "BF = {:.2f}".format(burstiness_factor_1))

#     yticks = [0, 0.2,  0.4,  0.6,  0.8, 1]
#     xticks = [0, 0.05,  0.1,  0.15,  0.2]


#     for ax in burst_axes:
#         ax.axvline(burst_threshold_scaled, color=axis_grey)
#         ax.set_ylim([0, 1])
#         ax.set_xlim([0, .23])
#         ax.set_yticks(yticks)
#         ax.set_xticks(xticks)
#         ax.set_ylabel("Frequency")
#         ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")

#     ax6.set_xlabel("Event duration (s)")


#     ax7.plot(scaled_g_BKs*1000, burstiness_factors_g_BK, marker=".")
#     ax7.set_xlabel(r"$\rightarrow g_{BK}$ (nS)")
#     ax7.set_ylabel("Burstiness")
#     ax7.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")
#     ax7.set_yticks(yticks)
#     ax7.set_ylim([-0.05, 1.05])
#     # ax7.set_xlim([0, scale_conductance(1)*1000])
#     # ax7.set_xlabel(r"$g_{BK}$ (mS/cm$^2$)")


#     ax8.plot(scaled_tau_BK, burstiness_factors_tau_BK, marker=".")
#     ax8.set_xlabel(r"$\tau_{BK}$ (ms)")
#     ax8.set_ylabel("Burstiness")
#     ax8.set_yticks(yticks)
#     ax8.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")
#     ax8.set_ylim([-0.05, 1.05])
#     ax8.set_xlim([2, 10])

#     ax7.text(label_x, label_y, r"\textbf{D}", transform=ax7.transAxes, fontsize=titlesize)
#     ax8.text(label_x, label_y, r"\textbf{E}", transform=ax8.transAxes, fontsize=titlesize)

#     plt.tight_layout()

#     plt.savefig("figure_1" + figure_format)



def figure_2():
    """
    Recreate figure 2 in Tabak et. al. 2011. Figure is saved as figure_1.png

    http://www.jneurosci.org/content/31/46/16855/tab-article-info
    """

    # g_bk => 0
    g_BK = scale_conductance(0)
    bins_0, binned_burstiness_factors_0 = robustness(g_BK=g_BK)

    # g => 0.5
    g_BK = scale_conductance(0.5)
    bins_05, binned_burstiness_factors_05 = robustness(g_BK=g_BK)

    # g_bk => 1
    g_BK = scale_conductance(1)
    bins_1, binned_burstiness_factors_1 = robustness(g_BK=g_BK)


    # Plotting
    plt.rcParams.update(params)

    fig, axes = plt.subplots(nrows=3, figsize=(figure_width, 1.5*figure_width))

    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]

    ax1.text(label_x, label_y, r"\textbf{A}", transform=ax1.transAxes, fontsize=titlesize)
    ax2.text(label_x, label_y, r"\textbf{B}", transform=ax2.transAxes, fontsize=titlesize)
    ax3.text(label_x, label_y, r"\textbf{C}", transform=ax3.transAxes, fontsize=titlesize)

    ax1.bar(bins_0[:-1], binned_burstiness_factors_0, width=(bins_0[1] - bins_0[0]), align="edge")
    title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(0)*1000) + r"$ (mS/cm$^2$) $\rightarrow 0$ nS"
    ax1.set_title(title, fontsize=titlesize)

    ax2.bar(bins_05[:-1], binned_burstiness_factors_05, width=(bins_05[1] - bins_05[0]), align="edge")
    title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(0.5)*1000) + r"$ (mS/cm$^2$) $\rightarrow 0.5$ nS"
    ax2.set_title(title, fontsize=titlesize)


    ax3.bar(bins_1[:-1], binned_burstiness_factors_1, width=(bins_1[1] - bins_1[0]), align="edge")
    title = r"$g_{BK} = " + "{:.2f}".format(scale_conductance(1)*1000) + r"$ (mS/cm$^2$) $\rightarrow 1$ nS"
    ax3.set_title(title, fontsize=titlesize)

    xticks = np.arange(0, 1.1, 0.2)

    for ax in [ax1, ax2, ax3]:
        ax.set_ylim([0, 350])
        ax.set_xticks(xticks)
        ax.set_ylabel("Number of models")
        ax.set_xlabel("Burstiness")
        ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")


    plt.tight_layout()

    plt.savefig("figure_2" + figure_format)



if __name__ == "__main__":
    figure_1()
    # figure_2()

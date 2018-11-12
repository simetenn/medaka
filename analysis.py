from __future__ import absolute_import, division, print_function, unicode_literals

import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib.patches import Patch

from medaka import rat, medaka_1, medaka_2
from burstiness import burstiness, duration, burst_threshold


# Simulation parameters
discard = 10000                # ms
simulation_time = 60000        # ms
noise_amplitude = 0.004        # nA
simulation_time_plot = 15000   # ms

# Plotting parameters
figure_width = 7.08
titlesize = 12
labelsize = 10
fontsize = 8
figure_format = ".eps"
label_x = -0.08
label_y = 1.08
axis_grey = (0.6, 0.6, 0.6)

figure_folder = "figures"
output_file = "spikes.txt"


# Set the random seed to increase reproducability
np.random.seed(10)


# Set default options for plotting
params = {
    "text.usetex": True,
    "text.latex.preamble": "\\usepackage{amsmath}, \\usepackage{amssymb}",
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


def calculate_frequency_bf(model, **parameters):
    """
    Calculate the frequency of a binned burstiness factor from a series of
    model evaluations.

    Parameters
    ----------
    model : func
        The model function.
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
    time, voltage = model(noise_amplitude=noise_amplitude,
                          discard=discard,
                          simulation_time=simulation_time,
                          **parameters)

    event_durations = duration(time, voltage)

    binned_durations, bins = np.histogram(event_durations,
                                          bins=nr_bins,
                                          range=hist_range)

    frequency = binned_durations/binned_durations.sum()

    burstiness_factor = burstiness(event_durations)

    return bins, frequency, burstiness_factor



def change_g_BK(model, original_g_BK, **parameters):
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
    scale_g_BKs = np.arange(0, 1.01, 0.025)

    burstiness_factors = []

    for scale_g_BK in scale_g_BKs:
        bins, frequency, burstiness_factor = \
            calculate_frequency_bf(model=model,
                                   g_BK=original_g_BK*scale_g_BK,
                                   **parameters)

        burstiness_factors.append(burstiness_factor)

    return scale_g_BKs, burstiness_factors




def comparison_noise():
    """
    Create the figure that compares RAT, MEDAKA 1 and MEDAKA 2 with noise,
    as well as the histogram of the frequency event durations.
    Figure saved as comparison_noise.
    """
    g_BK_scale = 0

    # Medaka 2
    time_0_medaka_2, V_0_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=4*3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time)

    bins_0_medaka_2, frequency_0_medaka_2, burstiness_factor_0_medaka_2 \
        = calculate_frequency_bf(medaka_2, g_BK=4*3.2e-4*g_BK_scale)

    # Medaka 1
    time_0_medaka_1, V_0_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time)

    bins_0_medaka_1, frequency_0_medaka_1, burstiness_factor_0_medaka_1 \
        = calculate_frequency_bf(medaka_1, g_BK=3.2e-4*g_BK_scale)

    # rat
    time_0_rat, V_0_rat = rat(noise_amplitude=noise_amplitude,
                              discard=discard,
                              g_BK=3.2e-4*g_BK_scale,
                              simulation_time=simulation_time)

    bins_0_rat, frequency_0_rat, burstiness_factor_0_rat \
        = calculate_frequency_bf(rat, g_BK=3.2e-4*g_BK_scale)


    g_BK_scale = 0.5

    # Medaka 2
    time_05_medaka_2, V_05_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                               discard=discard,
                                               g_BK=4*3.2e-4*g_BK_scale,
                                               simulation_time=simulation_time)

    bins_05, frequency_05, burstiness_factor_05 \
        = calculate_frequency_bf(medaka_2, g_BK=4*3.2e-4*g_BK_scale)

    # Medaka 1
    time_05_medaka_1, V_05_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                               discard=discard,
                                               g_BK=3.2e-4*g_BK_scale,
                                               simulation_time=simulation_time)

    bins_05_medaka_1, frequency_05_medaka_1, burstiness_factor_05_medaka_1 \
        = calculate_frequency_bf(medaka_1, g_BK=3.2e-4*g_BK_scale)


    # rat
    time_05_rat, V_05_rat = rat(noise_amplitude=noise_amplitude,
                                discard=discard,
                                g_BK=3.2e-4*g_BK_scale,
                                simulation_time=simulation_time)

    bins_05_rat, frequency_05_rat, burstiness_factor_05_rat \
        = calculate_frequency_bf(rat, g_BK=3.2e-4*g_BK_scale)


    g_BK_scale = 1

    # Medaka 2
    time_1_medaka_2, V_1_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=4*3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time)

    bins_1_medaka_2, frequency_1_medaka_2, burstiness_factor_1_medaka_2 = \
        calculate_frequency_bf(medaka_2, g_BK=4*3.2e-4*g_BK_scale,)

    # Medaka 1
    time_1_medaka_1, V_1_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time)

    bins_1_medaka_1, frequency_1_medaka_1, burstiness_factor_1_medaka_1 \
        = calculate_frequency_bf(medaka_1, g_BK=3.2e-4*g_BK_scale)

    # rat
    time_1_rat, V_1_rat = rat(noise_amplitude=noise_amplitude,
                              discard=discard,
                              g_BK=3.2e-4*g_BK_scale,
                              simulation_time=simulation_time)

    bins_1_rat, frequency_1_rat, burstiness_factor_1_rat \
        = calculate_frequency_bf(rat, g_BK=3.2e-4*g_BK_scale)

    # Calculate results for figure 1D

    # Medaka 2
    scaled_g_BKs_medaka_2, burstiness_factors_g_BK_medaka_2 \
        = change_g_BK(medaka_2, original_g_BK=4*3.2e-4)

    # Medaka 1
    scaled_g_BKs_medaka_1, burstiness_factors_g_BK_medaka_1 \
        = change_g_BK(medaka_1, original_g_BK=3.2e-4)

    # rat
    scaled_g_BKs_rat, burstiness_factors_g_BK_rat \
        = change_g_BK(rat, original_g_BK=3.2e-4)


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

    ax7 = plt.subplot(gs[3, 1:-1])


    voltage_axes = [ax1, ax2, ax3]
    burst_axes = [ax4, ax5, ax6]

    time_1_rat -= discard
    time_1_medaka_1 -= discard
    time_1_medaka_2 -= discard

    time_05_rat -= discard
    time_05_medaka_1 -= discard
    time_05_medaka_2 -= discard

    time_0_rat -= discard
    time_0_medaka_1 -= discard
    time_0_medaka_2 -= discard

    ax1.plot(time_1_rat, V_1_rat, color="tab:gray")
    ax1.plot(time_1_medaka_1, V_1_medaka_1, color="tab:blue")
    ax1.plot(time_1_medaka_2, V_1_medaka_2, color="tab:red")
    title = r"$1\cdot g_{\mathrm{BK}}$"
    ax1.set_title(title)
    ax1.text(label_x, label_y, r"\textbf{A}", transform=ax1.transAxes, fontsize=titlesize)


    ax2.plot(time_05_rat, V_05_rat, color="tab:gray")
    ax2.plot(time_05_medaka_1, V_05_medaka_1, color="tab:blue")
    ax2.plot(time_05_medaka_2, V_05_medaka_2, color="tab:red")
    title = r"$0.5\cdot g_{\mathrm{BK}}$"
    ax2.set_title(title)
    ax2.text(label_x, label_y, r"\textbf{B}", transform=ax2.transAxes, fontsize=titlesize)



    ax3.plot(time_0_rat, V_0_rat, color="tab:gray")
    ax3.plot(time_0_medaka_1, V_0_medaka_1, color="tab:blue")
    ax3.plot(time_0_medaka_2, V_0_medaka_2, color="tab:red")
    title = r"$0\cdot g_{\mathrm{BK}}$"
    ax3.set_title(title)
    ax3.text(label_x, label_y, r"\textbf{C}", transform=ax3.transAxes, fontsize=titlesize)
    ax3.set_xlabel("Time (ms)", fontsize=labelsize)

    for ax in voltage_axes:
        ax.set_ylabel("V (mV)")
        ax.set_ylim([-75, 40])
        ax.set_xlim([0, simulation_time_plot - discard])
        ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")



    ax4.bar(bins_1_medaka_1[:-1],
            frequency_1_medaka_1,
            width=(bins_1_medaka_1[1] - bins_1_medaka_1[0]),
            align="edge",
            color="tab:blue")
    ax4.text(100, 0.75, "BF = {:.2f}".format(burstiness_factor_1_medaka_1), color="tab:blue")

    ax4.bar(bins_1_medaka_2[:-1],
            frequency_1_medaka_2,
            width=(bins_1_medaka_2[1] - bins_1_medaka_2[0]),
            align="edge",
            color="tab:red")
    ax4.text(100, 0.6, "BF = {:.2f}".format(burstiness_factor_1_medaka_2), color="tab:red")

    ax4.bar(bins_1_rat[:-1],
            frequency_1_rat,
            width=(bins_1_rat[1] - bins_1_rat[0]),
            align="edge",
            color="tab:grey")
    ax4.text(100, 0.9, "BF = {:.2f}".format(burstiness_factor_1_rat), color="tab:grey")




    ax5.bar(bins_05_medaka_1[:-1],
            frequency_05_medaka_1,
            width=(bins_05_medaka_1[1] - bins_05_medaka_1[0]),
            align="edge",
            color="tab:blue")
    ax5.text(100, 0.75, "BF = {:.2f}".format(burstiness_factor_05_medaka_1), color="tab:blue")

    ax5.bar(bins_05[:-1],
            frequency_05,
            width=(bins_05[1] - bins_05[0]),
            align="edge",
            color="tab:red")
    ax5.text(100, 0.6, "BF = {:.2f}".format(burstiness_factor_05), color="tab:red")

    ax5.bar(bins_05_rat[:-1],
            frequency_05_rat,
            width=(bins_05_rat[1] - bins_05_rat[0]),
            align="edge",
            color="tab:grey")
    ax5.text(100, 0.9, "BF = {:.2f}".format(burstiness_factor_05_rat), color="tab:grey")


    ax6.bar(bins_0_rat[:-1],
            frequency_0_rat,
            width=(bins_0_rat[1] - bins_0_rat[0]),
            align="edge",
            color="tab:grey")
    ax6.text(100, 0.9, "BF = {:.2f}".format(burstiness_factor_0_rat), color="tab:grey")

    ax6.bar(bins_0_medaka_1[:-1],
            frequency_0_medaka_1,
            width=(bins_0_medaka_1[1] - bins_0_medaka_1[0]),
            align="edge",
            color="tab:blue")
    ax6.text(100, 0.75, "BF = {:.2f}".format(burstiness_factor_0_medaka_1), color="tab:blue")

    ax6.bar(bins_0_medaka_2[:-1],
            frequency_0_medaka_2,
            width=(bins_0_medaka_2[1] - bins_0_medaka_2[0]),
            align="edge",
            color="tab:red")
    ax6.text(100, 0.6, "BF = {:.2f}".format(burstiness_factor_0_medaka_2), color="tab:red")

    yticks = [0, 0.2,  0.4,  0.6,  0.8, 1]
    xticks = [0, 50,  100,  150,  200, 250]


    for ax in burst_axes:
        ax.axvline(burst_threshold, color=axis_grey)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, .3])
        ax.set_xticks(xticks)
        ax.set_ylabel("Frequency")
        ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")

    ax6.set_xlabel("Event duration (ms)")


    ax7.plot(scaled_g_BKs_rat, burstiness_factors_g_BK_rat, marker=".", color="tab:grey")
    ax7.plot(scaled_g_BKs_medaka_1, burstiness_factors_g_BK_medaka_1, marker=".", color="tab:blue")
    ax7.plot(scaled_g_BKs_medaka_2, burstiness_factors_g_BK_medaka_2, marker=".", color="tab:red")
    ax7.set_xlabel(r"$X \cdot g_{BK}$")
    ax7.set_ylabel("Burstiness")
    ax7.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")
    ax7.set_yticks(yticks)
    ax7.set_ylim([-0.05, 1.05])

    ax7.text(label_x, label_y, r"\textbf{D}", transform=ax7.transAxes, fontsize=titlesize)


    rat_legend = Patch(facecolor="tab:grey", label="RAT")
    medaka_1_legend = Patch(facecolor="tab:blue", label="MEDAKA 1")
    medaka_2_legend = Patch(facecolor="tab:red", label="MEDAKA 2")
    plt.legend(handles=[rat_legend, medaka_1_legend, medaka_2_legend],
               bbox_to_anchor=(0.8, 1),
               bbox_transform=plt.gcf().transFigure,
               ncol=3)


    plt.tight_layout()

    plt.subplots_adjust(top=0.91)

    plt.savefig(os.path.join(figure_folder, "comparison_noise" + figure_format))



def comparison_no_noise():
    """
    Create the figure that compares RAT, MEDAKA 1 and MEDAKA 2 without noise,
    with examples of the action potential shape.
    Figure saved as comparison_no_noise.
    """
    # g_bk => 0
    g_BK_scale = 0
    noise_amplitude = 0

    # Medaka 2
    time_0_medaka_2, V_0_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=4*3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time_plot)

    # Medaka 1
    time_0_medaka_1, V_0_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time_plot)

    # rat
    time_0_rat, V_0_rat = rat(noise_amplitude=noise_amplitude,
                              discard=discard,
                              g_BK=3.2e-4*g_BK_scale,
                              simulation_time=simulation_time_plot)


    # g_bk => 0
    g_BK_scale = 0.16

    # Medaka 2
    time_016_medaka_2, V_016_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                                 discard=discard,
                                                 g_BK=4*3.2e-4*g_BK_scale,
                                                 simulation_time=simulation_time_plot)

    # Medaka 1
    time_016_medaka_1, V_016_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                                 discard=discard,
                                                 g_BK=3.2e-4*g_BK_scale,
                                                 simulation_time=simulation_time_plot)

    # rat
    time_016_rat, V_016_rat = rat(noise_amplitude=noise_amplitude,
                                  discard=discard,
                                  g_BK=3.2e-4*g_BK_scale,
                                  simulation_time=simulation_time_plot)


    # g_bk => 0
    g_BK_scale = 0.2

    # Medaka 2
    time_02_medaka_2, V_02_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                               discard=discard,
                                               g_BK=4*3.2e-4*g_BK_scale,
                                               simulation_time=simulation_time_plot)

    # Medaka 1
    time_02_medaka_1, V_02_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                               discard=discard,
                                               g_BK=3.2e-4*g_BK_scale,
                                               simulation_time=simulation_time_plot)

    # rat
    time_02_rat, V_02_rat = rat(noise_amplitude=noise_amplitude,
                                discard=discard,
                                g_BK=3.2e-4*g_BK_scale,
                                simulation_time=simulation_time_plot)



    # g_bk => 0.5
    g_BK_scale = 0.5

    # Medaka 2
    time_05_medaka_2, V_05_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                               discard=discard,
                                               g_BK=4*3.2e-4*g_BK_scale,
                                               simulation_time=simulation_time_plot)

        # Medaka 1
    time_05_medaka_1, V_05_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                               discard=discard,
                                               g_BK=3.2e-4*g_BK_scale,
                                               simulation_time=simulation_time_plot)


    # rat
    time_05_rat, V_05_rat = rat(noise_amplitude=noise_amplitude,
                                discard=discard,
                                g_BK=3.2e-4*g_BK_scale,
                                simulation_time=simulation_time_plot)

    # g_bk => 1
    g_BK_scale = 1

    # Medaka 2
    time_1_medaka_2, V_1_medaka_2 = medaka_2(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=4*3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time_plot)


    # Medaka 1
    time_1_medaka_1, V_1_medaka_1 = medaka_1(noise_amplitude=noise_amplitude,
                                             discard=discard,
                                             g_BK=3.2e-4*g_BK_scale,
                                             simulation_time=simulation_time_plot)

    # rat
    time_1_rat, V_1_rat = rat(noise_amplitude=noise_amplitude,
                              discard=discard,
                              g_BK=3.2e-4*g_BK_scale,
                              simulation_time=simulation_time_plot)


    # Plotting
    plt.rcParams.update(params)

    fig = plt.figure(figsize=(figure_width, 1.35*figure_width))
    gs = plt.GridSpec(5, 6)
    ax1 = plt.subplot(gs[0, :-2])
    ax2 = plt.subplot(gs[1, :-2])
    ax3 = plt.subplot(gs[2, :-2])
    ax4 = plt.subplot(gs[3, :-2])
    ax5 = plt.subplot(gs[4, :-2])

    ax6 = plt.subplot(gs[0, -2:])
    ax7 = plt.subplot(gs[1, -2:])
    ax8 = plt.subplot(gs[2, -2:])
    ax9 = plt.subplot(gs[3, -2:])
    ax10 = plt.subplot(gs[4, -2:])


    voltage_axes = [ax1, ax2, ax3, ax4, ax5]
    selection_axes = [ax6, ax7, ax8, ax9, ax10]

    time_1_rat -= discard
    time_1_medaka_1 -= discard
    time_1_medaka_2 -= discard

    time_05_rat -= discard
    time_05_medaka_1 -= discard
    time_05_medaka_2 -= discard

    time_02_rat -= discard
    time_02_medaka_1 -= discard
    time_02_medaka_2 -= discard

    time_016_rat -= discard
    time_016_medaka_1 -= discard
    time_016_medaka_2 -= discard

    time_0_rat -= discard
    time_0_medaka_1 -= discard
    time_0_medaka_2 -= discard

    ax1.plot(time_1_rat, V_1_rat, color="tab:gray")
    ax1.plot(time_1_medaka_1, V_1_medaka_1, color="tab:blue")
    ax1.plot(time_1_medaka_2, V_1_medaka_2, color="tab:red")
    title = r"$1\cdot g_{\mathrm{BK}}$"
    ax1.set_title(title)
    ax1.text(label_x, label_y, r"\textbf{A1}", transform=ax1.transAxes, fontsize=titlesize)


    ax2.plot(time_05_rat, V_05_rat, color="tab:gray")
    ax2.plot(time_05_medaka_1, V_05_medaka_1, color="tab:blue")
    ax2.plot(time_05_medaka_2, V_05_medaka_2, color="tab:red")
    title = r"$0.5\cdot g_{\mathrm{BK}}$"
    ax2.set_title(title)
    ax2.text(label_x, label_y, r"\textbf{B1}", transform=ax2.transAxes, fontsize=titlesize)


    ax3.plot(time_02_rat, V_02_rat, color="tab:gray")
    ax3.plot(time_02_medaka_1, V_02_medaka_1, color="tab:blue")
    ax3.plot(time_02_medaka_2, V_02_medaka_2, color="tab:red")
    title = r"$0.2\cdot g_{\mathrm{BK}}$"
    ax3.set_title(title)
    ax3.text(label_x, label_y, r"\textbf{C1}", transform=ax3.transAxes, fontsize=titlesize)


    ax4.plot(time_016_rat, V_016_rat, color="tab:gray")
    ax4.plot(time_016_medaka_1, V_016_medaka_1, color="tab:blue")
    ax4.plot(time_016_medaka_2, V_016_medaka_2, color="tab:red")
    title = r"$0.16\cdot g_{\mathrm{BK}}$"
    ax4.set_title(title)
    ax4.text(label_x, label_y, r"\textbf{D1}", transform=ax4.transAxes, fontsize=titlesize)


    ax5.plot(time_0_rat, V_0_rat, color="tab:gray")
    ax5.plot(time_0_medaka_1, V_0_medaka_1, color="tab:blue")
    ax5.plot(time_0_medaka_2, V_0_medaka_2, color="tab:red")
    title = r"$0\cdot g_{\mathrm{BK}}$"
    ax5.set_title(title)
    ax5.text(label_x, label_y, r"\textbf{E1}", transform=ax5.transAxes, fontsize=titlesize)
    ax5.set_xlabel("Time (ms)", fontsize=labelsize)


    def get_index(time, start):
        delta_time = 500

        return np.where((time > start) & (time < start + delta_time))[0]


    label_x_shift = -0.1

    # Selections
    with open(output_file, "w") as output:
        # 1*g_bk
        index = get_index(time_1_rat, 2500)
        tmp_time = time_1_rat[index] - time_1_rat[index].min()
        ax6.plot(tmp_time, V_1_rat[index], color="tab:gray")
        output.write("1*g_bk RAT: " +  str(V_1_rat[index].max()) + "\n")

        index = get_index(time_1_medaka_1, 2400)
        tmp_time = time_1_medaka_1[index] - time_1_medaka_1[index].min()
        ax6.plot(tmp_time, V_1_medaka_1[index], color="tab:blue")
        output.write("1*g_bk MEDAKA 1: " + str(V_1_medaka_1[index].max()) + "\n")

        index = get_index(time_1_medaka_2, 2320)
        tmp_time = time_1_medaka_2[index] - time_1_medaka_2[index].min()
        ax6.plot(tmp_time, V_1_medaka_2[index], color="tab:red")
        ax6.text(label_x + label_x_shift, label_y, r"\textbf{A2}", transform=ax6.transAxes, fontsize=titlesize)
        ax6.set_title("Action potential shape")
        output.write("1*g_bk MEDAKA 2: " + str(V_1_medaka_2[index].max()) + "\n\n")


        # 0.5*g_bk
        index = get_index(time_05_rat, 2450)
        tmp_time = time_05_rat[index] - time_05_rat[index].min()
        ax7.plot(tmp_time, V_05_rat[index], color="tab:gray")
        output.write("0.5*g_bk RAT: " + str(V_05_rat[index].max()) + "\n")

        index = get_index(time_05_medaka_1, 2570)
        tmp_time = time_05_medaka_1[index] - time_05_medaka_1[index].min()
        ax7.plot(tmp_time, V_05_medaka_1[index], color="tab:blue")
        output.write("0.5*g_bk MEDAKA 1: " + str(V_05_medaka_1[index].max()) + "\n")

        index = get_index(time_05_medaka_2, 2000)
        tmp_time = time_05_medaka_2[index] - time_05_medaka_2[index].min()
        ax7.plot(tmp_time, V_05_medaka_2[index], color="tab:red")
        ax7.text(label_x + label_x_shift, label_y, r"\textbf{B2}", transform=ax7.transAxes, fontsize=titlesize)
        output.write("0.5*g_bk MEDAKA 2: " + str(V_05_medaka_2[index].max()) + "\n\n")


        # 0.2*g_bk
        index = get_index(time_02_rat, 900)
        tmp_time = time_02_rat[index] - time_02_rat[index].min()
        ax8.plot(tmp_time, V_02_rat[index], color="tab:gray")
        output.write("0.2*g_bk RAT: " + str(V_02_rat[index].max()) + "\n")

        index = get_index(time_02_medaka_1, 1050)
        tmp_time = time_02_medaka_1[index] - time_02_medaka_1[index].min()
        ax8.plot(tmp_time, V_02_medaka_1[index], color="tab:blue")
        output.write("0.2*g_bk MEDAKA 1: " + str(V_02_medaka_1[index].max()) + "\n")

        index = get_index(time_02_medaka_2, 1280)
        tmp_time = time_02_medaka_2[index] - time_02_medaka_2[index].min()
        ax8.plot(tmp_time, V_02_medaka_2[index], color="tab:red")
        ax8.text(label_x + label_x_shift, label_y, r"\textbf{C2}", transform=ax8.transAxes, fontsize=titlesize)
        output.write("0.2*g_bk MEDAKA 2: " + str(V_02_medaka_2[index].max()) + "\n\n")


        # 0.16*g_bk
        index = get_index(time_016_rat, 3050)
        tmp_time = time_016_rat[index] - time_016_rat[index].min()
        ax9.plot(tmp_time, V_016_rat[index], color="tab:gray")
        output.write("0.16*g_bk RAT: " + str(V_016_rat[index].max()) + "\n")

        index = get_index(time_016_medaka_1, 3350)
        tmp_time = time_016_medaka_1[index] - time_016_medaka_1[index].min()
        ax9.plot(tmp_time, V_016_medaka_1[index], color="tab:blue")
        output.write("0.16*g_bk MEDAKA 1: " + str(V_016_medaka_1[index].max()) + "\n")

        index = get_index(time_016_medaka_2, 3100)
        tmp_time = time_016_medaka_2[index] - time_016_medaka_2[index].min()
        ax9.plot(tmp_time, V_016_medaka_2[index], color="tab:red")
        ax9.text(label_x + label_x_shift, label_y, r"\textbf{D2}", transform=ax9.transAxes, fontsize=titlesize)
        output.write("0.16*g_bk MEDAKA 2: " + str(V_016_medaka_2[index].max()) + "\n\n")


        # 0*g_bk
        index = get_index(time_0_rat, 2000)
        tmp_time = time_0_rat[index] - time_0_rat[index].min()
        ax10.plot(tmp_time, V_0_rat[index], color="tab:gray")
        output.write("0*g_bk RAT: " + str(V_0_rat[index].max()) + "\n")

        index = get_index(time_0_medaka_1, 2200)
        tmp_time = time_0_medaka_1[index] - time_0_medaka_1[index].min()
        ax10.plot(tmp_time, V_0_medaka_1[index], color="tab:blue")
        output.write("0*g_bk MEDAKA 1: " + str(V_0_medaka_1[index].max()) + "\n")

        index = get_index(time_0_medaka_2, 3850)
        tmp_time = time_0_medaka_2[index] - time_0_medaka_2[index].min()
        ax10.plot(tmp_time, V_0_medaka_2[index], color="tab:red")
        output.write("0*g_bk MEDAKA 2:" + str(V_0_medaka_2[index].max()) + "\n")


    ax10.set_xlabel("Relative time (ms)", fontsize=labelsize)
    ax10.text(label_x + label_x_shift, label_y, r"\textbf{E2}", transform=ax10.transAxes, fontsize=titlesize)


    for ax in voltage_axes:
        ax.set_ylabel("V (mV)")
        ax.set_ylim([-75, 40])
        ax.set_xlim([0, simulation_time_plot - discard])
        ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")


    for ax in selection_axes:
        ax.set_ylim([-75, 40])
        ax.set_xlim([0, 300])
        ax.tick_params(axis="both", which="major", labelsize=fontsize, labelcolor="black")



    rat_legend = Patch(facecolor="tab:grey", label="RAT")
    medaka_1_legend = Patch(facecolor="tab:blue", label="MEDAKA 1")
    medaka_2_legend = Patch(facecolor="tab:red", label="MEDAKA 2")
    plt.legend(handles=[rat_legend, medaka_1_legend, medaka_2_legend],
               bbox_to_anchor=(0.8, 1),
               bbox_transform=plt.gcf().transFigure,
               ncol=3)


    plt.tight_layout()

    plt.subplots_adjust(top=0.91, wspace=1)

    plt.savefig(os.path.join(figure_folder, "comparison_no_noise" + figure_format))



if __name__ == "__main__":
    comparison_noise()
    comparison_no_noise()
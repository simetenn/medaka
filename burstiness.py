from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import uncertainpy as un


onset_threshold = 0.55          # fraction of the normalized voltage
burst_threshold = 60            # in ms
end_threshold = -0.1            # Relative to the onset_threshold
min_spike_amplitude = 10        # mV


def duration(time, voltage):
    """
    Calculate the duration of a series of events in a voltage trace.

    Parameters
    ----------
    time : array_like
        Time array, in ms.
    voltage : array_like
        Voltage array from a neuron.

    Returns
    -------
    duration : numpy.array
        A numpy array containing the duration of each spike (in ms) in `voltage`.

    Notes
    -----
    Normalizes the voltage before the spikes are calculated.
    """
    # Find spikes in the normalized voltage trace
    spikes = un.features.Spikes(time,
                                voltage,
                                threshold=onset_threshold,
                                end_threshold=end_threshold,
                                trim=False,
                                normalize=True,
                                min_amplitude=min_spike_amplitude)

    # Calculate the duration of each spike
    duration = []
    for spike in spikes:
        duration.append(spike.time[-1] - spike.time[0])

    return np.array(duration)



def burstiness(durations, burst_threshold=burst_threshold):
    """
    Calculate the burstiness from a series of event durations by finding the
    fraction of events that is above `burst_threshold`.

    Parameters
    ----------
    durations : array
        A numpy array containing the duration of each spike in `voltage`.
    burst_threshold : float, optional
        Default is 60 ms.

    Returns
    -------
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """

    durations = np.array(durations)

    if not durations.size:
        return None

    # Find all events greater than the burst_threshold
    burst_occurrences = durations > burst_threshold

    # Calculate the fraction of these events
    burstiness_factor = np.sum(burst_occurrences)/len(burst_occurrences)

    return burstiness_factor



def bursting(time, spikes, info):
    """
    If the model has bursts or not. Is one if the model has at least one spike
    with duration longer than the burstiness threshold, else it is zero.

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
    found_bursts : int
        Is one if the model has at least one spike with duration longer than
        the burstiness threshold, else it is zero.
    """
    found_bursts = 0
    for spike in spikes:
        if (spike.time[-1] - spike.time[0]) > burst_threshold:
            found_bursts = 1
            break

    return None, found_bursts


def spiking(time, spikes, info):
    """
    If the model has spikes or not. Is one if the model has at least one spike
    with duration shorter than the burstiness threshold, else it is zero.

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
        Is one if the model has at least one spike with duration shorter than
        the burstiness threshold, else it is zero.
    """
    found_spike = 0
    for spike in spikes:
        if (spike.time[-1] - spike.time[0]) <= burst_threshold:
            found_spike = 1
            break

    return None, found_spike


def APs(time, spikes, info):
    """
    If the model have action potentials or not. Is one if the model has at
    least one spike, else it is zero.

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
    found_APs : int
        Is one if the model has at least one spike, else it is zero.
    """
    found_APs = 0
    if spikes.nr_spikes > 0:
        found_APs = 1

    return None, found_APs
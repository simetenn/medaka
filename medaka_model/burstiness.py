from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import uncertainpy as un


onset_threshold = 0.55          # fraction of the normalized voltage
burst_threshold = 60            # in ms
end_threshold = -0.1            # Relative to the onset_threshold

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
                                normalize=True)

    # Calculate the duration of each spike
    duration = []
    for spike in spikes:
        duration.append(spike.time[-1] - spike.time[0])

    return np.array(duration)



def burstiness(durations, burst_threshold=60):
    """
    Calculate the burstiness from a series of event durations by finding the
    fraction of events that is above `burst_threshold`.

    Parameters
    ----------
    durations : array
        A numpy array containing the duration of each spike in `voltage`.
    burst_threshold : float, optional
        Default is 70 ms.

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



def burstiness_factor(time, spikes, info):
    """
    Calculate the burstiness factor from a uncertainpy.Spikes object.
    Compatible with uncertainpy.SpikingFeatures.

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
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """
    # Calculate the duration of each spike
    durations = []
    for spike in spikes:
        durations.append(spike.time[-1] - spike.time[0])

    burstiness_factor = burstiness(durations, burst_threshold=burst_threshold)

    return None, burstiness_factor



def bursting(time, spikes, info):
    """
    Compatible with uncertainpy.SpikingFeatures.

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
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """
    for spike in spikes:
        if (spike.time[-1] - spike.time[0]) > burst_threshold:
            return None, 1

    return None, 0


def spiking(time, spikes, info):
    """
    Compatible with uncertainpy.SpikingFeatures.

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
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """
    for spike in spikes:
        if (spike.time[-1] - spike.time[0]) < burst_threshold:
            return None, 1

    return None, 0


def APs(time, spikes, info):
    """
    Compatible with uncertainpy.SpikingFeatures.

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
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """
    APs = 0
    if spikes.nr_spikes > 0:
        APs = 1

    return None, APs



def bursts_spikes_aps(time, spikes, info):
    """
    Compatible with uncertainpy.SpikingFeatures.

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
    burstiness_factor : float, None
        The fraction of events that is above `burst_threshold`. Returns None
        if there are no events.
    """
    voltage = spikes.V
    time = spikes.time

    durations = duration(time, voltage)

    if not durations.size:
        bursts = -1
    elif np.any(durations > burst_threshold):
        bursts = 1
    else:
        bursts = 0


    return None, bursts

# def burstiness_efel(time, voltage, info):
#     """
#     Calculate the burstiness factor from a uncertainpy.Spikes object.
#     Compatible with uncertainpy.EfelFeatures.

#     Parameters
#     ----------
#     time : array_like
#         Time array, in ms.
#     voltage : array_like
#         Voltage array from a neuron.
#     info : dictionary
#         Not used, but is required as input

#     Returns
#     -------
#     time : None
#     burstiness_factor : float, None
#         The fraction of events that is above `burst_threshold`. Returns None
#         if there are no events.
#     """
#     event_durations = duration(time, voltage)

#     burstiness_factor = burstiness(event_durations, burst_threshold=burst_threshold)

#     return None, burstiness_factor
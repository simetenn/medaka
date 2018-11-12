from __future__ import absolute_import, division, print_function, unicode_literals


import numpy as np

import neuron
nrn = neuron.h

# time step
dt = 0.25 # ms

A = 3.1415927e-6 # cm^2

c_NEURON = 10*1e-6/A
alpha_NEURON = 0.0015*A*10**6


g_l = 6.37e-5              # S/cm^2
g_Ca_rat = 0               # S/cm^2
g_SK = 6.37e-4             # S/cm^2
g_Na = 0.07                # S/cm^2
g_Ca = 2e-4                # S/cm^2
g_K = 9.55e-4              # S/cm^2
tau_n = 30
v_f = -20

def create_soma(g_l=g_l,
                g_K=g_K,
                g_SK=g_SK,
                g_BK=0,
                g_Na=g_Na,
                g_Ca=g_Ca,
                g_Ca_rat=g_Ca_rat,
                E_leak=-45,
                tau_BK=5,
                tau_n=tau_n,
                v_f=v_f):
    """
    Create the soma of a neuron.

    Parameters
    ----------
    g_l : float, optional
        The leak conductance, in S/cm^2. Default is 6.37e-5 S/cm^2.
    g_K : float, optional
        The maximal conductance of K channels, in S/cm^2. Default is
        9.55e-4 S/cm^2.
    g_SK : float, optional
        The maximal conductance of SK channels, in S/cm^2. Default is
        6.37e-4 S/cm^2.
    g_BK : float, optional
        The maximal conductance of BK channels, in S/cm^2. Default is 0 S/cm^2.
    g_Na: float, optional
        The maximal conductance of Na channels, in S/cm^2. Default is 0.07 S/cm^2.
    g_Ca : float, optional
        The maximal conductance of Ca channels (medaka version), in S/cm^2.
        Default is 2e-4 S/cm^2.
    g_Ca_rat : float, optional
        The maximal conductance of Ca channels (tabak version), in S/cm^2. Default is
        0 S/cm^2.
    E_leak : float, optional
        Reversal potential for the leak current, in mV. Default is -45 mV.
    tau_BK : float, optional
        Time constant of the BK channel, in ms. Default is 5 ms.
    tau_n : float, optional
        Time constant of n (activation of I_K), in ms. Default is 30 ms.
    v_f : float, optional
        Voltage value at the midpoint of f (activation of I_BK), in mV. Default
        is -20 mV.

    Returns
    -------
    soma : soma NEURON object
        The soma NEURON object.
    """
    nrn('forall delete_section()')
    soma = nrn.Section('soma')
    soma.L = 10                            # um; stored as a float number
    soma.diam = 10                         # um
    soma.nseg = 1                          # stored as an integer

    for sec in nrn.allsec():
        sec.Ra = 100
        sec.cm = c_NEURON

        sec.insert('pas')
        sec.insert("kdrt")                 # From Tabak 2011
        sec.insert("Cadt")                 # From Tabak 2011:
        sec.insert("ihvat")                # From Halnes 2011:
        sec.insert("sk")                   # From Halnes 2011
        sec.insert("bk")                   # From Tabak 2011:

        sec.insert("naxm")                 # From Migle, but updated to Medaka data
        sec.insert("ihva")                 # From Halnes 2011, but updated to Medaka data

        sec.ena = 50                       # Reversal potential for sodium
        sec.ek = -75                       # Reversal potential for potassium
        sec.eCa = 60                       # Reversal potential for calcium


        for seg in sec:
            seg.ftau_bk = tau_BK
            seg.alpha_Cadt = alpha_NEURON
            seg.g_pas = g_l
            seg.e_pas = E_leak
            seg.gkdrbar_kdrt = g_K # med kun kdrt & nax er 0.001 og 0.05 fine tall
            seg.ghvat_ihvat = g_Ca_rat
            seg.gskbar_sk = g_SK
            seg.gbk_bk = g_BK
            seg.taun_kdrt = tau_n
            seg.vf_bk = v_f

            # TODO Check these
            seg.gbar_naxm = g_Na
            seg.pcabar_ihva = g_Ca

    return soma


def insert_current_clamp(input_site, duration=5000, delay=0, amplitude=0):
    """
    Inserts a current clamp in the neuron model.

    Parameters
    ----------
    input_site : neuron.Segment
        Where to place the current clamp. Example: soma(0.5), where 0.5 means 'center',
        0 would mean start, and 1 would mean at the end of the segment in question.
    duration : {float, int}, optional
        Duration of stimulus in ms. Default is 5000 ms.
    delay: float, optional
        Delay of stimulus in ms. Default is 0 ms.
    amplitude:
        Amplitude of stimulus in nA. Default is 0 nA.

    Returns
    -------
    stim : NEURON object current clamp
        The NEURON object current clamp. This must be returned, otherwise it is
        lost.
    """
    stim = nrn.IClamp(input_site)
    stim.delay = delay
    stim.dur = duration
    stim.amp = amplitude

    return stim


def run_simulation(soma, simulation_time=5000, noise_amplitude=0):
    """
    Runs the NEURON simulation.

    Parameters
    ----------
    record_site : neuron.Segment
        Where to record membrane potential from. Example: soma(0.5), where 0.5
        means 'center', 0 would mean start, and 1 would mean at the end of the
        segment in question.
    simulation_time : {float, int}, optional
        Simulation time in ms. Default is 5000 ms.
    noise_amplitude : float, optional
        The amplitude of the noise added to the model, in nA. If 0, no noise is added.
        Note that the model uses adaptive timesteps if there is no noise,
        and fixed timesteps with dt = 0.25 if there is noise. Default is 0.

    Returns
    -------
    time : array
        Time array for the simulation.
    voltage : array
        Voltage array for the simulation.
    """
    rec_t, rec_v = record(soma(0.5))

    cvode = nrn.CVode()

    if noise_amplitude == 0:
        cvode.active(1)

        nrn.finitialize(-60)
        neuron.init()

        neuron.run(simulation_time)

    else:
        cvode.active(0)

        noise_stim = insert_current_clamp(soma(0.5), duration=simulation_time)

        nrn.dt = dt
        nrn.finitialize(-60)
        neuron.init()

        n_steps = int(np.ceil(simulation_time/nrn.dt)) + 1
        noise = noise_amplitude*np.random.normal(size=n_steps)/np.sqrt(nrn.dt)

        # Add noise
        i = 0
        while nrn.t < simulation_time:
            noise_stim.amp = noise[i]
            nrn.fadvance()

            i += 1

    return np.array(rec_t), np.array(rec_v)


def record(record_site):
    """
    Set up time and voltage recordings.

    Parameters
    ----------
    record_site : neuron.Segment
        Where to record voltage.

    Returns
    -------
    rec_t : neuron.hocObject : A Neuron Vector object.
        A Neuron Vector object for the time.
    rec_v : neuron.hocObject : A Neuron Vector object.
        A Neuron Vector object for the voltage.
    """
    rec_t = nrn.Vector()
    rec_t.record(nrn._ref_t)

    rec_v = nrn.Vector()
    rec_v.record(record_site._ref_v)

    return rec_t, rec_v


def run_general_model(g_l=g_l,
                      g_K=g_K,
                      g_SK=g_SK,
                      g_BK=0,
                      g_Na=g_Na,
                      g_Ca=g_Ca,
                      g_Ca_rat=g_Ca_rat,
                      E_leak=-45,
                      tau_BK=5,
                      tau_n=tau_n,
                      v_f=v_f,
                      simulation_time=5000,
                      noise_amplitude=0,
                      stimulus_amplitude=0,
                      discard=0):
    """
    Run the model of pituitary cells. Default set of parameters is the same
    as for medaka_1.

    Parameters
    ----------
    g_l : float, optional
        The leak conductance, in S/cm^2. Default is 6.37e-5 S/cm^2.
    g_K : float, optional
        The maximal conductance of K channels, in S/cm^2. Default is
        9.55e-4 S/cm^2.
    g_SK : float, optional
        The maximal conductance of SK channels, in S/cm^2. Default is
        6.37e-4 S/cm^2.
    g_BK : float, optional
        The maximal conductance of BK channels, in S/cm^2. Default is 0 S/cm^2.
    g_Na: float, optional
        The maximal conductance of Na channels, in S/cm^2. Default is 0.07 S/cm^2.
    g_Ca : float, optional
        The maximal conductance of Ca channels (medaka version), in S/cm^2.
        Default is 2e-4 S/cm^2.
    g_Ca_rat : float, optional
        The maximal conductance of Ca channels (tabak version), in S/cm^2. Default is
        0 S/cm^2.
    E_leak : float, optional
        Reversal potential for the leak current, in mV. Default is -45 mV.
    tau_BK : float, optional
        Time constant of the BK channel, in ms. Default is 5 ms.
    tau_n : float, optional
        Time constant of n (activation of I_K), in ms. Default is 30 ms.
    v_f : float, optional
        Voltage value at the midpoint of f (activation of I_BK), in mV. Default
        is -20 mV.
    discard : {float, int}, optional
        The first ms of the simulation to be discarded. Default is 0 ms.
    simulation_time : {float, int}, optional
        Simulation time in ms. Default is 5000 ms.
    noise_amplitude : float, optional
        The amplitude of the noise added to the model, in nA. If 0, no noise is
        added. Note that the model uses adaptive timesteps if there is no noise,
        and fixed timesteps with dt = 0.25 if there is noise. Default is 0.
    stimulus_amplitude : float, optional
        The amplitude of the stimulus added to the model, in nA. Default is 0.

    Returns
    -------
    time : array
        Time array for the simulation.
    voltage : array
        Voltage array for the simulation.
    """
    soma = create_soma(g_l=g_l,
                       g_K=g_K,
                       g_SK=g_SK,
                       g_BK=g_BK,
                       g_Na=g_Na,
                       g_Ca=g_Ca,
                       g_Ca_rat=g_Ca_rat,
                       E_leak=E_leak,
                       tau_BK=tau_BK,
                       tau_n=tau_n,
                       v_f=v_f)

    stim = insert_current_clamp(soma(0.5),
                                duration=simulation_time,
                                amplitude=stimulus_amplitude)

    time, voltage = run_simulation(soma,
                                   simulation_time=simulation_time,
                                   noise_amplitude=noise_amplitude)

    return time[time > discard], voltage[time > discard]



def rat(g_l=g_l,
        g_K=g_K,
        g_SK=g_SK,
        g_BK=0,
        g_Ca=6.37e-4,
        E_leak=-50,
        tau_BK=5,
        tau_n=tau_n,
        v_f=v_f,
        simulation_time=5000,
        noise_amplitude=0,
        stimulus_amplitude=0,
        discard=0):
    """
    Neuron model of pituitary cells in rat, reimplementation of Tabak et. al. 2011:

    http://www.jneurosci.org/content/31/46/16855/tab-article-info

    g_l : float, optional
        The leak conductance, in S/cm^2. Default is 6.37e-5 S/cm^2.
    g_K : float, optional
        The maximal conductance of K channels, in S/cm^2. Default is
        9.55e-4 S/cm^2.
    g_SK : float, optional
        The maximal conductance of SK channels, in S/cm^2. Default is
        6.37e-4 S/cm^2.
    g_BK : float, optional
        The maximal conductance of BK channels, in S/cm^2. Default is 0 S/cm^2.
    g_Na: float, optional
        The maximal conductance of Na channels, in S/cm^2. Default is 0.07 S/cm^2.
    g_Ca : float, optional
        The maximal conductance of Ca channels (tabak version), in S/cm^2. Default is
        6.37e-4 S/cm^2.
    E_leak : float, optional
        Reversal potential for the leak current, in mV. Default is -50 mV.
    tau_BK : float, optional
        Time constant of the BK channel, in ms. Default is 5 ms.
    tau_n : float, optional
        Time constant of n (activation of I_K), in ms. Default is 30 ms.
    v_f : float, optional
        Voltage value at the midpoint of f (activation of I_BK), in mV. Default
        is -20 mV.
    discard : {float, int}, optional
        The first ms of the simulation to be discarded. Default is 0 ms.
    simulation_time : {float, int}, optional
        Simulation time in ms. Default is 5000 ms.
    noise_amplitude : float, optional
        The amplitude of the noise added to the model, in nA. If 0, no noise is
        added. Note that the model uses adaptive timesteps if there is no noise,
        and fixed timesteps with dt = 0.25 if there is noise. Default is 0.
    stimulus_amplitude : float, optional
        The amplitude of the stimulus added to the model, in nA. Default is 0.


    Returns
    -------
    time : array
        Time array for the simulation.
    voltage : array
        Voltage array for the simulation.
    """
    time, voltage = run_general_model(g_l=g_l,
                                      g_K=g_K,
                                      g_SK=g_SK,
                                      g_BK=g_BK,
                                      g_Na=0,
                                      g_Ca=0,
                                      g_Ca_rat=g_Ca,
                                      E_leak=E_leak,
                                      tau_BK=tau_BK,
                                      tau_n=tau_n,
                                      v_f=v_f,
                                      simulation_time=simulation_time,
                                      noise_amplitude=noise_amplitude,
                                      stimulus_amplitude=stimulus_amplitude,
                                      discard=discard)

    return time, voltage



def medaka_1(g_l=g_l,
             g_K=g_K,
             g_SK=g_SK,
             g_BK=0,
             g_Na=g_Na,
             g_Ca=g_Ca,
             E_leak=-45,
             tau_BK=5,
             tau_n=tau_n,
             v_f=v_f,
             simulation_time=5000,
             noise_amplitude=0,
             stimulus_amplitude=0,
             discard=0):
    """
    Run the model of pituitary cells. Default set of parameters is the same
    as for medaka_1.

    Parameters
    ----------
    g_l : float, optional
        The leak conductance, in S/cm^2. Default is 6.37e-5 S/cm^2.
    g_K : float, optional
        The maximal conductance of K channels, in S/cm^2. Default is
        9.55e-4 S/cm^2.
    g_SK : float, optional
        The maximal conductance of SK channels, in S/cm^2. Default is
        6.37e-4 S/cm^2.
    g_BK : float, optional
        The maximal conductance of BK channels, in S/cm^2. Default is 0 S/cm^2.
    g_Na: float, optional
        The maximal conductance of Na channels, in S/cm^2. Default is 0.07 S/cm^2.
    g_Ca : float, optional
        The maximal conductance of Ca channels (medaka version), in S/cm^2.
        Default is 2e-4 S/cm^2.
    E_leak : float, optional
        Reversal potential for the leak current, in mV. Default is -45 mV.
    tau_BK : float, optional
        Time constant of the BK channel, in ms. Default is 5 ms.
    tau_n : float, optional
        Time constant of n (activation of I_K), in ms. Default is 30 ms.
    v_f : float, optional
        Voltage value at the midpoint of f (activation of I_BK), in mV. Default
        is -20 mV.
    discard : {float, int}, optional
        The first ms of the simulation to be discarded. Default is 0 ms.
    simulation_time : {float, int}, optional
        Simulation time in ms. Default is 5000 ms.
    noise_amplitude : float, optional
        The amplitude of the noise added to the model, in nA. If 0, no noise is
        added. Note that the model uses adaptive timesteps if there is no noise,
        and fixed timesteps with dt = 0.25 if there is noise. Default is 0.
    stimulus_amplitude : float, optional
        The amplitude of the stimulus added to the model, in nA. Default is 0.

    Returns
    -------
    time : array
        Time array for the simulation.
    voltage : array
        Voltage array for the simulation.
    """
    time, voltage = run_general_model(g_l=g_l,
                                      g_K=g_K,
                                      g_SK=g_SK,
                                      g_BK=g_BK,
                                      g_Na=g_Na,
                                      g_Ca=g_Ca,
                                      g_Ca_rat=0,
                                      E_leak=E_leak,
                                      tau_BK=tau_BK,
                                      tau_n=tau_n,
                                      v_f=v_f,
                                      simulation_time=simulation_time,
                                      noise_amplitude=noise_amplitude,
                                      stimulus_amplitude=stimulus_amplitude,
                                      discard=discard)
    return time, voltage


def medaka_2(g_l=g_l,
             g_K=g_K*1.4,
             g_SK=g_SK*3,
             g_BK=4*3.2e-4,
             g_Na=g_Na,
             g_Ca=g_Ca,
             E_leak=-45,
             tau_BK=5,
             tau_n=5,
             v_f=-15,
             simulation_time=5000,
             noise_amplitude=0,
             stimulus_amplitude=0,
             discard=0):
    """
    Medaka 2 neuron model of medaka cells in fish. Tuned to experimental results.

    Parameters
    ----------
    g_l : float, optional
        The leak conductance, in S/cm^2. Default is 6.37e-5 S/cm^2.
    g_K : float, optional
        The maximal conductance of K channels, in S/cm^2. Default is
        1.337e-3 S/cm^2.
    g_SK : float, optional
        The maximal conductance of SK channels, in S/cm^2. Default is
        1.92e-3 S/cm^2.
    g_BK : float, optional
        The maximal conductance of BK channels, in S/cm^2. Default is 1.28e-3 S/cm^2.
    g_Na: float, optional
        The maximal conductance of Na channels, in S/cm^2. Default is 0.07 S/cm^2.
    g_Ca : float, optional
        The maximal conductance of Ca channels (medaka version), in S/cm^2.
        Default is 2e-4 S/cm^2.
    E_leak : float, optional
        Reversal potential for the leak current, in mV. Default is -45 mV.
    tau_BK : float, optional
        Time constant of the BK channel, in ms. Default is 5 ms.
    tau_n : float, optional
        Time constant of n (activation of I_K), in ms. Default is 5 ms.
    v_f : float, optional
        Voltage value at the midpoint of f (activation of I_BK), in mV. Default
        is -15 mV.
    discard : {float, int}, optional
        The first ms of the simulation to be discarded. Default is 0 ms.
    simulation_time : {float, int}, optional
        Simulation time in ms. Default is 5000 ms.
    noise_amplitude : float, optional
        The amplitude of the noise added to the model, in nA. If 0, no noise is
        added. Note that the model uses adaptive timesteps if there is no noise,
        and fixed timesteps with dt = 0.25 if there is noise. Default is 0.
    stimulus_amplitude : float, optional
        The amplitude of the stimulus added to the model, in nA. Default is 0.

    Returns
    -------
    time : array
        Time array for the simulation.
    voltage : array
        Voltage array for the simulation.

    Notes
    -----
    Compared to Medaka 1, this model has changed the parameters g_K, g_BK, tau_n
    and v_f. The kinetics of the models are unchanged.
    """
    time, voltage = run_general_model(g_l=g_l,
                                      g_K=g_K,
                                      g_SK=g_SK,
                                      g_BK=g_BK,
                                      g_Na=g_Na,
                                      g_Ca=g_Ca,
                                      g_Ca_rat=0,
                                      E_leak=E_leak,
                                      tau_BK=tau_BK,
                                      tau_n=tau_n,
                                      v_f=v_f,
                                      simulation_time=simulation_time,
                                      noise_amplitude=noise_amplitude,
                                      stimulus_amplitude=stimulus_amplitude,
                                      discard=discard)

    return time, voltage



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    time, V = medaka_1()

    plt.style.use("seaborn-darkgrid")
    plt.figure()
    plt.plot(time, V)
    plt.savefig("voltage.png")


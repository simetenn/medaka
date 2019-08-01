#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 09:28:40 2018

@author: geih
"""

# Initialize the model and defining default options
import numpy as np
import neuron
nrn = neuron.h
import uncertainpy as un
import pylab as plt

plt.interactive(1)
plt.show()


def create_soma(g_l=0.064e-3,
                e_pas=-45,
                g_K=1.3e-3,
                pcabar_ihva=0.2e-3,
                g_SK=0.32e-3,
                g_BK=1e-3,
                gbar_naxm=0.07,
                tau_BK=3,
                tau_K=5):

    nrn('forall delete_section()')
    soma = nrn.Section('soma')
    soma.L = 10                            # um; stored as a float number
    soma.diam = 10                         # um
    soma.nseg = 1                          # stored as an integer

    for sec in nrn.allsec():
        sec.Ra = 100
        sec.cm = 1.0
        sec.insert('pas')
        sec.insert("kdrt")                 # From Tabak 2011
        sec.insert("Cadt")                 # From Tabak 2011:
        sec.insert("sk")                   # From Halnes 2011
        sec.insert("bk")                   # From Tabak 2011:
        sec.insert("naxm")                 # From Migle, but updated to Medaka data
        sec.insert("ihva")                 # From Halnes 2011, but updated to Medaka data

        sec.ena = 50                       # Reversal potential for sodium
        sec.ek = -75                       # Reversal potential for potassium
        sec.eCa = 60                       # Reversal potential for calcium


        for seg in sec:
            seg.alpha_Cadt = 1.51e-2
            seg.g_pas = g_l
            seg.e_pas = e_pas
            seg.gkdrbar_kdrt = g_K # med kun kdrt & nax er 0.001 og 0.05 fine tall
            seg.pcabar_ihva = pcabar_ihva
            seg.gskbar_sk = g_SK
            seg.gbk_bk = g_BK
            seg.gbar_naxm = gbar_naxm
            seg.ntau_bk = tau_BK
            seg.taun_kdrt = tau_K
            seg.AA_bk = 1.21

    return soma



def insert_current_clamp(input_site, simulation_time=5000):
    stim = nrn.IClamp(input_site)
    stim.delay = 0
    stim.dur = simulation_time
    stim.amp = 0

    return stim



def record(record_site):
    rec_t = nrn.Vector()
    rec_t.record(nrn._ref_t)
    rec_v = nrn.Vector()
    rec_v.record(record_site._ref_v)
    rec_ca = nrn.Vector()
    rec_ca.record(record_site._ref_Cai)

    return rec_t, rec_v, rec_ca


def run_simulation(record_site, stim, simulation_time=5000, noise_amplitude=0):
    rec_t, rec_v, rec_ca = record(record_site)
    cvode = nrn.CVode()

    if noise_amplitude == 0:
        cvode.active(1)

        nrn.finitialize(-60)
        neuron.init()

        neuron.run(simulation_time)

    else:
        cvode.active(0)

        nrn.dt = 0.25
        nrn.finitialize(-60)
        neuron.init()

        n_steps = int(np.ceil(simulation_time/nrn.dt)) + 1
        noise = noise_amplitude*np.random.normal(size=n_steps)/np.sqrt(nrn.dt)

        # Add noise
        i = 0
        while nrn.t < simulation_time:
            stim.amp = noise[i]
            nrn.fadvance()

            i += 1

    return np.array(rec_t), np.array(rec_v), np.array(rec_ca)



simulation_time = 16000        # in ms
noise_amplitude = 0*0.001            # in mV

epasval = -45 # Tabak value -50
glval = 2e-5 # Tabak value 0.64e-4
Kval = 4.18e-4 # Tabak value 9.55e-4
SKval	= 4e-4 # Tabak-value 6.4e-4: increase to reduce firing rate
tauKval = 5 #
CaFishval = 6.25e-5 # Tabak-value was 6.4e-4, but different kinetics
Naval = 2.19e-2 # Not part of Tabak-model
BKdef = 3.13e-4 #
tauBKval = 3 # Tau for BK-near. Duncan had 20 ms


def medaka(e_pas=epasval,
           g_l=glval,
           g_K=Kval,
           g_Ca=CaFishval,
           g_SK=SKval,
           g_BK=BKdef,
           g_Na=Naval,
           tau_BK=tauBKval,
           tau_K=tauKval,
           simulation_time=simulation_time,
           noise_amplitude=noise_amplitude):

    soma = create_soma(e_pas=e_pas,
                       g_l=g_l,
                       g_K=g_K,
                       pcabar_ihva=g_Ca,
                       g_SK=g_SK,
                       g_BK=g_BK,
                       gbar_naxm=g_Na,
                       tau_BK=tau_BK,
                       tau_K=tau_K)

    stim = insert_current_clamp(soma(0.5), simulation_time=simulation_time)

    t, v, ca = run_simulation(soma(0.5),
                              stim,
                              simulation_time=simulation_time,
                              noise_amplitude=noise_amplitude)

    v = v[(t>2850) & (t<12850)]
    t = t[(t>2850) & (t<12850)]
    t = t-t[0]

    return t, v



if __name__ == "__main__":
    for ii in range(0, 7):
        if ii == 0:
            BKval = BKdef # BK near
        if ii == 1:
            BKval = 0.25*BKdef # Tabak-value up to 3.2e-4
        if ii == 2:
            BKval = 0.16*BKdef # Tabak-value up to 3.2e-4
        if ii == 3:
            BKval = 0.15*BKdef # Tabak-value up to 3.2e-4
        if ii == 4:
            BKval = 0.13*BKdef #0.152*BKdef # Tabak-value up to 3.2e-4
        if ii == 5:
            BKval = 0
        if ii == 6:
            BKval = BKdef # Tabak-value up to 3.2e-4
            Naval = 0 # TTX effect

        soma = create_soma(e_pas=epasval,
                g_l=glval,
                g_K=Kval,
                pcabar_ihva=CaFishval,
                g_SK=SKval,
                g_BK=BKval,
                gbar_naxm=Naval,
                tau_BK = tauBKval,
                tau_K = tauKval)


        stim = insert_current_clamp(soma(0.5), simulation_time=simulation_time)

        t, v, ca = run_simulation(soma(0.5), stim,
                        simulation_time=simulation_time,
                        noise_amplitude=noise_amplitude)

        v = v[(t>2850) & (t<12850)]
        ca = ca[(t>2850) & (t<12850)]
        ca = ca*1e6 # convert to nM
        t = t[(t>2850) & (t<12850)]
        t = t-t[0]

        if ii == 0:
            t0, v0, ca0 = t, v, ca
        if ii == 1:
            t1, v1, ca1 = t, v, ca
        if ii == 2:
            t2, v2, ca2 = t, v, ca
        if ii == 3:
            t3, v3, ca3 = t, v, ca
        if ii == 4:
            t4, v4, ca4 = t, v, ca
        if ii == 5:
            t5, v5, ca5 = t, v, ca
        if ii == 6:
            t6, v6, ca6 = t, v, ca

    fig = plt.figure(23)
    # plt.rc('text', usetex=True)
    #plt.gcf().text(0.02, 0.9, "Fish Gon.", fontsize=14)


    ax1 = fig.add_subplot(6,3,(1,2), xlabel="", ylabel="$V_m$ [mV]", title = "A1")
    ax2 = fig.add_subplot(6,3,3, xlabel="", ylabel="", title = "A2")

    ax3 = fig.add_subplot(6,3,(4,5), xlabel="", ylabel="$V_m$ [mV]", title = "B1")
    ax4 = fig.add_subplot(6,3,6, xlabel="", ylabel="", title = "B2")

    ax5 = fig.add_subplot(6,3,(7,8), xlabel="", ylabel="$V_m$ [mV]", title = "C1")
    ax6 = fig.add_subplot(6,3,9, xlabel="", ylabel="", title = "C2")

    ax7 = fig.add_subplot(6,3,(10,11), xlabel="", ylabel="$V_m$ [mV]", title = "D1")
    ax8 = fig.add_subplot(6,3,12, xlabel="", ylabel="", title = "D2")

    ax9 = fig.add_subplot(6,3,(13,14), xlabel="", ylabel="$V_m$ [mV]", title = "E1")
    ax10 = fig.add_subplot(6,3,15, xlabel="", ylabel="", title = "E2")

    ax11 = fig.add_subplot(6,3,(16,17), xlabel="$t$ [ms]", ylabel="$V_m$ [mV]", title = "F1")
    ax12 = fig.add_subplot(6,3,18, xlabel="$t$ [ms]", ylabel="", title = "F2")



    ax1.plot(t0, v0, 'tab:blue',label='control')

    ax1.plot(t6, v6, 'tab:red',label='$g_{Na}=0$')

    #ax2.plot(t0, v0, 'tab:blue', label='control')
    ax2.plot(t0, v0, 'tab:blue', label='control')
    ax2.plot(t6, v6, 'tab:red', label='$g_{Na}=0$')
    ax1.legend(loc = 1)


    ax3.plot(t1, v1, 'tab:blue', label='$g_{BK}\cdot0.25$')

    ax4.plot(t1, v1, 'tab:blue', label='$g_{BK}\cdot$ 0.25')
    ax3.legend(loc = 1)

    ax5.plot(t2, v2, 'tab:blue', label='$g_{BK}\cdot$ 0.16')

    ax6.plot(t2, v2, 'tab:blue', label='$g_{BK}\cdot$ 0.16')
    ax5.legend(loc = 1)

    ax7.plot(t3, v3, 'tab:blue', label='$g_{BK}\cdot$ 0.15')

    ax8.plot(t3, v3, 'tab:blue', label='$g_{BK}\cdot$ 0.15')
    ax7.legend(loc = 1)

    ax9.plot(t4, v4, 'tab:blue', label='$g_{BK}\cdot$ 0.13')
    ax10.plot(t4, v4, 'tab:blue', label='$g_{BK}\cdot$ 0.13')
    ax9.legend(loc = 1)

    ax11.plot(t5, v5, 'tab:blue', label='$g_{BK}\cdot$ 0')
    ax12.plot(t5, v5, 'tab:blue', label='$g_{BK}\cdot$ 0')
    ax11.legend(loc = 1)

    #ax1.set_xlim([3000,13000])
    #ax3.set_xlim([3000,13000])
    #ax5.set_xlim([3000,13000])
    #ax7.set_xlim([3000,13000])
    #ax9.set_xlim([3000,13000])
    #ax11.set_xlim([3000,13000])

    ax2.set_xlim([0,150])
    ax4.set_xlim([0,150])
    ax6.set_xlim([0,150])
    ax8.set_xlim([0,150])
    ax10.set_xlim([0,150])
    ax12.set_xlim([0,150])


    #ax3.set_xlim([0,160])
    #ax4.set_xlim([0,400])

    ax1.set_ylim([-70,10])
    ax2.set_ylim([-70,10])
    ax3.set_ylim([-70,10])
    #ax4.set_ylim([-70,10])

    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xticks([])
    ax4.set_xticks([])
    ax5.set_xticks([])
    ax6.set_xticks([])
    ax7.set_xticks([])
    ax8.set_xticks([])
    ax9.set_xticks([])
    ax10.set_xticks([])



    fig.subplots_adjust(hspace=0.4)
    fig.subplots_adjust(wspace=0.4)

    #plt.legend(loc=4, frameon=False)
    plt.savefig('medakasim.png')

    fig = plt.figure(4)
    # plt.rc('text', usetex=True)
    plt.gcf().text(0.02, 0.9, "MEDAKA 2$", fontsize=14)

    ax1 = fig.add_subplot(1,1,1, xlabel="time [ms]", ylabel="Calcium [nM]", title = "A) Calcium")
    ax1.plot(t3, ca3, 'tab:grey', label = '$g_{BK} \cdot 0$')
    ax1.plot(t2, ca2, 'tab:green', label = '$g_{BK} \cdot 0.16$', linewidth = 3.0)
    ax1.plot(t1, ca1, 'y', label = '$g_{BK} \cdot 0.18$')
    ax1.plot(t4, ca4, 'tab:pink', label = '$g_{BK} \cdot 0.20$')
    ax1.plot(t0, ca0, 'tab:blue', label = '$g_{BK} \cdot 1$', linewidth = 3.0)
    ax1.plot(t5, ca5, 'tab:red', label = '$g_{Na} \cdot 0$')
    ax1.legend(loc=1)
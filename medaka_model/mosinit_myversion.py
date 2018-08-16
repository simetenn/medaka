#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 2018


@author: geih
"""

from __future__ import division
from scipy.stats import norm
from neuron import h, load_mechanisms
from numpy import trapz
cvode = h.CVode()
cvode.active(1)
#cvode.maxstep(0.2)
h.load_file('stdlib.hoc')

import numpy as np

import pylab as plt
plt.interactive(1) 


import neuron
nrn = neuron.h

import matplotlib.cm as cm

celsius = 20.0 # temperature

def return_soma():
    """
    Makes a cell model.
    :return: soma NEURON object, using single compartment
    """
    nrn('forall delete_section()')
    soma = nrn.Section('soma')
    soma.L = 10  # um; stored as a float number
    soma.diam = 10  # um
    soma.nseg = 1  # stored as an integer
    
    for sec in nrn.allsec():
        sec.insert('pas')
        sec.Ra = 100
        sec.cm = 1
        sec.insert("kdrt") # From Tabak 2011
        sec.insert("naxm") # From Migle, but fitted to medakadata
        sec.insert("Cad") # From Halnes 2011
        sec.insert("ihva") # From Halnes 2011, but fitted to medakadata
        sec.insert("sk") # From Halnes 2011, but adapted to model in Tabak 2011
        sec.insert("bk") # From Tabak 2011

        sec.ena = 50 # Reversal potential for sodium
        sec.ek = -65 # Reversal potential for potassium
        sec.eCa = 100 # Reversal potential for calcium

        for seg in sec:
            seg.g_pas = 0.00003
            seg.e_pas = -55
#            seg.gbar_naxm = 0.05

            seg.gkdrbar_kdrt = 0.001 # med kun kdrt & nax er 0.001 og 0.05 fine tall
            seg.gbar_naxm = 0.06
            seg.pcabar_ihva = 1e-5
            seg.gskbar_sk	= 2e-5
            seg.gbk_bk = 3e-4
    return soma

def insert_current_clamp(input_site, arrivetime):
    """
    Inserts a current clamp in the neuron model
    :param input_site: Where to place the current clamp. Example: soma(0.5), where 0.5 means 'center',
           0 would mean start, and 1 would mean at the end of the segment in question.
    :return: The NEURON object current clamp. This must be returned, otherwise it is lost.
    """
    stim = nrn.IClamp(input_site)
    stim.delay = arrivetime
#    stim.amp = 0.0007
    stim.amp = 0.0007
    stim.dur = 3000
    return stim


def run_simulation(record_site):
    """
    Runs the NEURON simulation
    :param record_site: Where to record membrane potential from. Example: soma(0.5), where 0.5 means 'center',
           0 would mean start, and 1 would mean at the end of the segment in question.
    :return: Time and voltage numpy arrays
    """
    rec_t = nrn.Vector()
    rec_t.record(nrn._ref_t)
    rec_v = nrn.Vector()
    rec_v.record(record_site._ref_v)
    rec_ca = nrn.Vector()
    rec_ca.record(record_site._ref_Cai)
    neuron.h.dt = 2**-3
    nrn.finitialize(-55)
    neuron.init()
    neuron.run(3000)
    return np.array(rec_t), np.array(rec_v), np.array(rec_ca)


def myoutput_1():

    mysoma = return_soma()
    stim = insert_current_clamp(mysoma(0.5),500)
#    stim2 = insert_current_clamp(mysoma(0.5),1500)

    myt, myv, myca = run_simulation(mysoma(0.5))

    fig = plt.figure(1)
    ax1 = fig.add_subplot(111, xlabel="Time [ms]", ylabel="Voltage [mV]")
    ax1.plot(myt, myv, 'b', label='control')
    plt.legend(loc=4, frameon=False)

    plt.savefig('medakasim.png')


    fig = plt.figure(2)
    ax1 = fig.add_subplot(111, xlabel="Time [ms]", ylabel="Calcium [mM]")
    ax1.plot(myt, myca, 'r', label='more L')
    plt.legend(loc=4, frameon=False)

    plt.savefig('medakasim.png')


if __name__ == '__main__':
    myoutput_1()

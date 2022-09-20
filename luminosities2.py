#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 12:40:11 2022

@author: trinitystenhouse

Here is the module for calculations and plotting. The functions defined here
are what you will be using to get an output from your simulations.
"""
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from Geometries2 import *
import math
from scipy.stats import poisson

#%%
def lum_thresh(X_ref, X_E, S_d, S_ref, L_ref):
    """
    calculates threshold luminosity of an event using the reference 
    values found in DUNE paper

    Parameters
    ----------
    X_ref : FLOAT FOR REF CROSS SECTION
    X_E : ARRAY OF CROSS SECTIONS
        neutrino cross section array of our event
    S_d : NUMBER/ARRAY
        emission surface geometry ie sphere or torus of our event
    S_ref : FLOAT FOR REF GEOMETRY
    L_ref : FLOAT FOR REF LUMINOSITY 
        ie ref luminosity at event not at earth

    Returns
    -------
    L_t : ARRAY OF THRESHOLD LUMINOSITIES
        threshold luminosity of event at each cross section

    """
    
    L_t = []
    for i in range(0, len(X_E)):
        L_ti = L_ref * np.array(( X_ref / X_E[i] ) * ( S_d / S_ref )) #find detector area 
        L_t.append(L_ti)
        
    return L_t

#%%Piecewise Fit

def piecewise_linear(x, x0, x1, b, k1, k2, k3):
    #plots piecewise fit 
    
    condlist = [x < x0, (x >= x0) & (x < x1), x >= x1]
    funclist = [lambda x: k1*x + b, lambda x: k1*x + b + k2*(x-x0), lambda x: k1*x + b + k2*(x-x0) + k3*(x - x1)]
    return np.piecewise(x, condlist, funclist)

#%%Threshold Mechanisms

def threshold(L, E, param):
    #calculates how many events are aboove the threshold curve, by comparing to
    #the value of the fit at each energy and luminosity
    
    events_count = 0
    
    for i in range(0, len(L)):
        f = piecewise_linear(E[i], *param)
        if f <= L[i]:
            events_count += 1
        else:
            pass
        
    return events_count

def thresholdlog(L, L_t):
    #calculates how many events are aboove the threshold curve, by comparing to
    #the value of the fit at each energy and luminosity
    events_count = 0
    
    for i in range(0, len(L)):
        if L_t[i] < L[i]:
            events_count += 1
        else:
            pass
    print(events_count)
    return events_count

#%% Fitting Mechanisms

def fitting1(refdata, eventdata, L_t, variables):
    #standard graph, no log axes
    name = variables[0]
    t = variables[1]
    xlim = variables[2]
    ylim = variables[3]
    p0 = variables[4]
    
    param , cov = curve_fit(piecewise_linear, refdata[0], L_t, p0)
    
    #plotting the graph with all data points and threshold curve
    fig, ax = plt.subplots(figsize=(8, 8))
    plot = plt.scatter(eventdata[0], eventdata[3], c = t, cmap=cm.get_cmap('spring'))
    plt.errorbar(refdata[0], L_t, label = 'threshold', fmt = 'x', color = 'turquoise')
    plt.plot(refdata[0], piecewise_linear(refdata[0], *param), color = 'turquoise')
    plt.legend()
    cb = plt.colorbar(plot, ax=ax, label = 'Time (ms)')
    #ax.set_xscale('log')
    plt.title(name)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Luminosity (ergs/s)')
    plt.show()
    
    return param

def fittinglog(refdata, eventdata, L_t, name):
    #single dataset log graph - both x and y axes are log
    #plotting the graph with all data points and threshold curve
    fig, ax = plt.subplots(figsize=(8, 8))
    plot = plt.scatter(eventdata[0], eventdata[3], color = 'magenta')
    plt.errorbar(refdata[0], L_t, label = 'threshold', fmt = 'x', color = 'turquoise')
    plt.plot(refdata[0], L_t, color = 'turquoise')
    plt.legend()
    ax.set_yscale('log')
    plt.semilogx()
    plt.title(name)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Luminosity (ergs/s)')
    plt.show()
    
def fitmulti(refdata, eventdata, L_t, variables):
    #standard graph, no log axes, 2 datasets plotted
    #plotting the graph with all data points and threshold curve
    name = variables[0]
    t1 = variables[1]
    xlim = variables[2]
    ylim = variables[3]
    p0 = variables[4]
    t2 = variables[5]
    
    param , cov = curve_fit(piecewise_linear, refdata[0], L_t, p0)
    yerr = 10e50
    
    fig, ax = plt.subplots(figsize=(14, 8))
    plot_1 = plt.scatter(eventdata[0], eventdata[3], c = t1, cmap=cm.get_cmap('spring'))
    plot_2 = plt.scatter(eventdata[4], eventdata[7], c = t2, cmap=cm.get_cmap('winter'))
    plt.errorbar(refdata[0], L_t, label = 'threshold', yerr = yerr, fmt = 'x', color = 'turquoise')
    plt.plot(refdata[0], piecewise_linear(refdata[0], *param), color = 'turquoise')
    plt.legend()
    cb = plt.colorbar(plot_1, ax=ax, label = 'Time Irrotational (ms)')
    cb2 = plt.colorbar(plot_2, ax=ax, label = 'Time Corotational (ms)')
    plt.title(name)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel('Energy (MeV)')
    plt.ylabel( 'Luminosity (ergs/s)')
    plt.show()
    
    return param

def fittinglogmulti(refdata, eventdata, L_t, name):
    #2 dataset log graph, both axes are log 
    #plotting the graph with all data points and threshold curve
    fig, ax = plt.subplots(figsize=(8, 8))
    plot = plt.scatter(eventdata[0], eventdata[3], color = 'magenta', label = '0rads inclination')
    plot = plt.scatter(eventdata[4], eventdata[7], color = 'mediumpurple', label = '1 degree off $\pi / 2$rads inclination')
    plt.errorbar(refdata[0], L_t, label = 'threshold', fmt = 'x', color = 'turquoise')
    plt.plot(refdata[0], L_t, color = 'turquoise')
    plt.legend()
    ax.set_yscale('log')
    plt.semilogx()
    plt.title(name)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Luminosity (ergs/s)')
    plt.show()
    
#%% Full Information Function 
#plots a graph of your choice for your data, and returns a proportion of events above threshold

def full(refdata, eventdata, variables, fittype):
    """
    Parameters
    ----------
    refdata : array
        reference dataset [energy, X_ref, S_ref, L_ref]
    eventdata : array
        event(s) dataset 
        [E_v, X_E, S_d, Earth_Lum] or [E_v1, X_E1, S_d1, Earth_Lum1, E_v2, X_E2, S_d2, Earth_Lum2]
    variables : array
        [name, time, xlim, ylim, p0] and add t2 if relevant 
        all fit types need at least object in position 0 (name)
    fittype : string
        type of fitting you want for your dataset
        mine were made specifically for my datasets, so it may be necessary to edit
            some parameters ie if you want time included in 'log' type

    Returns
    -------
    counter : float or array of floats
        proportion of events above threshold

    """
    
    name = variables[0]
    
    #find the threshold luminosities for our dataset
    L_t = lum_thresh(refdata[1], eventdata[1], eventdata[2], refdata[2], refdata[3])
    #threshold is min L you see a few neutrinos
    
    if fittype == 'piecewise':
        #plot graph
        param = fitting1(refdata, eventdata, L_t, variables)
        #calc number of events above threshold
        events_count = threshold(eventdata[3], eventdata[0], param)
        percentage_events = events_count / len(eventdata[3]) * 100
        print(percentage_events, '% events above threshold')
        counter = percentage_events/100
    elif fittype == 'multi':
        param = fitmulti(refdata, eventdata, L_t, variables)
        #calc number of events above threshold
        events_count_1 = threshold(eventdata[3], eventdata[0], param)
        events_count_2 = threshold(eventdata[7], eventdata[4], param)
        percentage_events_1 = events_count_1 / len(eventdata[3]) * 100
        percentage_events_2 = events_count_2 / len(eventdata[3]) * 100
        counter = [percentage_events_1/100, percentage_events_2/100]
        print(percentage_events_1, '% events above threshold for event 1', \
              percentage_events_2, '% events above threshold for event 2')
    elif fittype == 'log':
        fittinglog(refdata, eventdata, L_t, name)
        #calc number of events above threshold
        events_count = thresholdlog(eventdata[3], L_t)
        percentage_events = events_count / len(eventdata[3]) * 100
        print(percentage_events, '% events above threshold')
        counter = percentage_events/100
    elif fittype == 'multilog':
        fittinglogmulti(refdata, eventdata, L_t, name)
        #calc number of events above threshold
        events_count_1 = thresholdlog(eventdata[3], L_t)
        events_count_2 = thresholdlog(eventdata[7], L_t)
        percentage_events_1 = events_count_1 / len(eventdata[3]) * 100
        percentage_events_2 = events_count_2 / len(eventdata[7]) * 100
        print(percentage_events_1, '% events above threshold for event 1')
        print(percentage_events_2, '% events above threshold for event 2')
        counter = [percentage_events_1/100, percentage_events_2/100]
    else:
        print('Not sure how to fit that yet!')
        
    return counter





#%%

from scipy.stats import poisson
from scipy.stats import norm
import matplotlib.mlab as mlab

def rate_events(X_E, Earth_Lum):
    """
    calculates rate of events being observed
    from cross section and luminosity

    """
    rate_events = []
    for i in range(0, len(Earth_Lum)):
        rate_events.append(X_E[i] * Earth_Lum[i])
    return rate_events # * detector area, detector efficiency

def conversion_factor(X, L, t, no_events):
    #converts the rate to the actual number of events compared
    #to the reference luminosities and events counts
    
    rate = X * L
    ratext = rate * t
    cf = no_events / ratext # = 1 / rate    
    return cf


def v_energy(spectrum, no_neutrinos):
    """
    chooses a random energy value for the neutrino based on 
    imported energy spectrum of merger
    
    also smooths spectrum so there are no gaps

    """ 
    recon_spectrum = []
    
    for j in range(0, len(no_neutrinos)):
        for i in range(0, int(no_neutrinos[j] * 10000)): 
            n = ((random.randint(int(spectrum[j] * 10000 - 5000), (int(spectrum[j] * 10000 + 5000))))/ 10000)
            recon_spectrum.append(n)
        
    return recon_spectrum

def energy_hist(t, X_ref, X_E, L_ref, L_Earth, no_events_ref, events_count, \
                t_ref, spectrum, bwc, name):
    """
    plots a histogram made via monte carlo simulation with a 
    fitted gaussian overlay
    """
    if events_count == 0:
        print('No events will be seen')
        
    else:
        cf = conversion_factor(X_ref, L_ref, t_ref, no_events_ref)
        print(cf, 'cf')
        actual_events = np.array(rate_events(X_E, L_Earth)) * cf
        
        
        energy_list = np.array(v_energy(spectrum, actual_events))
            #assigns a random energy to every neutrino 
                        #that causes an event 
        
        #plot histogram of energies
        plt.hist(energy_list, density=True, bins = 50, alpha = 0.25, color = 'turquoise', label = 'Raw Data')  
        plt.ylabel('Events')
        plt.title(name +' Histogram')
        plt.xlabel('Energies (MeV)')
        plt.legend()
        plt.show()
        
        return energy_list
    

#%% determining references for scale

def events(refdata, eventdata, Earth_Lum_ref, X_actual, t_ref, no_events_ref, percent, name):
    
    cf = conversion_factor(refdata[1], Earth_Lum_ref, t_ref, no_events_ref)
    actual_events = np.array(rate_events(X_actual, eventdata[3])) * cf * percent

    #prop_list = proportion_per_bin(actual_events)
    
    plt.title('Histogram 1')
    hist, bins, idk = plt.hist(eventdata[0], bins = len(eventdata[0]), weights = actual_events)
    plt.show()
    
    
    histactual =  []
    for i in range(0, len(eventdata[0])):
        histactual.append(int(hist[i]) * actual_events[i])
    
    if 10e29 <= np.mean(histactual):
        plt.title('Histogram for Poisson Fluctuations for '  
                  + name)
        hist, bins, idk = plt.hist(eventdata[0], weights = actual_events, bins = 10, color = 'turquoise', alpha = 0.25, label = 'Events')
        plt.xlabel('Energy (MeV)')
        plt.legend()
        plt.ylabel('Events')
        plt.show()
    else:
        s2 = np.random.poisson(lam=(actual_events), size = len(eventdata[0]))
        #poisson approximates to gaussian for high event numbers so this is fine to do
        
        plt.title('Histogram for Poisson Fluctuations for '  
                  + name)
        hist2, bins, idk = plt.hist(eventdata[0], weights = (s2), bins = 10, color = 'turquoise', alpha = 0.25, label = 'Events')
        plt.xlabel('Energy (MeV)')
        plt.legend()
        plt.ylabel('Events')
        plt.show()


#%% Graph of Distance against Proportion of events above threshold
def distgraph(refdata, eventdata, area, Earth_Lum, name, var, fittype):
    """

    Parameters
    ----------
    refdata : array
        reference [energy, cross section, area, luminosity]
    eventdata : array
        event's [energy, cross section, area, luminosity at Earth]
    area : float
        new area created from input distance array
    Earth_Lum : float
        earth luminosity calculated using input distance from distance array
    name : string
        graph name
    var : array
        any other variables needed to calculate fit
        variables in order [name, time, xlim, ylim, p0]
    fittype : string
        specifies which fit method to use

    Returns
    -------
    counter : float or array of floats
        proportion of events above threshold

    """
    
    #find the threshold luminosities for our dataset
    L_t = lum_thresh(refdata[1], eventdata[1], eventdata[2], refdata[2], refdata[3])
    #threshold is min L you see a few neutrinos
    
    eventdata = [eventdata[0], eventdata[1], eventdata[2], Earth_Lum]
    
    if fittype == 'piecewise':
        #plot graph
        p0=var[4]
        param , cov = curve_fit(piecewise_linear, refdata[0], L_t, p0)
        #calc number of events above threshold
        events_count = threshold(Earth_Lum, eventdata[0], param)
        percentage_events = events_count / len(Earth_Lum) * 100
        print(percentage_events, '% events above threshold')
        counter = percentage_events/100
    elif fittype == 'log':
        #calc number of events above threshold
        events_count = thresholdlog(Earth_Lum, L_t)
        percentage_events = events_count / len(Earth_Lum) * 100
        print(percentage_events, '% events above threshold')
        counter = percentage_events/100
    else:
        print('Not sure how to fit that yet!')

    return counter


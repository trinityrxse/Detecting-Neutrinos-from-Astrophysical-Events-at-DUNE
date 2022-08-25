#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 12:40:11 2022

@author: trinitystenhouse
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
    X_E : ARRAY OF CROSS SECTIONS
        neutrino cross section array of our event
    S_d : NUMBER/ARRAY
        emission surface geometry ie sphere or torus of our event

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

#%%


def piecewise_linear(x, x0, x1, b, k1, k2, k3):
    #plots piecewise fit 
    
    condlist = [x < x0, (x >= x0) & (x < x1), x >= x1]
    funclist = [lambda x: k1*x + b, lambda x: k1*x + b + k2*(x-x0), lambda x: k1*x + b + k2*(x-x0) + k3*(x - x1)]
    return np.piecewise(x, condlist, funclist)

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


def plot(energy, E_v, X_ref, X_E, S_d, S_ref, L_ref, Earth_Lum, name, t, xlim, ylim, p0):
    """
    

    Parameters
    ----------
    energy : array
        reference energy used to plot threshold curve
    E_v : array
        neutrino energy for the event in question
    X_ref : float
        reference cross section used to calculate threshold curve
    X_E : array
        cross section of neutrinos at the energies of the E_v array
    S_d : float
        surface area of event
    S_ref : float
        surface area of detector
    L_ref : float
        reference luminosity used to calculate threshold curve
    Earth_Lum : array
        neutrino luminosities of event in question at Earth 
        ie scaled depending on geometry
    name : string
        list element for title of graphs
    t : array
        times of events
    xlim : tuple
        domain of graph
    ylim : tuple
        range of graph 

    Returns
    -------
    events_count : float
        number of events above threshold

    """
    #find the threshold luminosities for our dataset
    L_t = lum_thresh(X_ref, X_E, S_d, S_ref, L_ref)
    #threshold is min L you see a few neutrinos
    
    
    p0 = p0
    param , cov = curve_fit(piecewise_linear, energy, L_t, p0)
    #yerr = 10e50
    
    #calc number of events above threshold
    events_count = threshold(Earth_Lum, E_v, param)

    #plotting the graph with all data points and threshold curve
    fig, ax = plt.subplots(figsize=(8, 8))
    plot = plt.scatter(E_v, Earth_Lum, c = t, cmap=cm.get_cmap('spring'))
    plt.errorbar(energy, L_t, label = 'threshold', fmt = 'x', color = 'turquoise')
    plt.plot(energy, piecewise_linear(energy, *param), color = 'turquoise')
    plt.legend()
    cb = plt.colorbar(plot, ax=ax, label = 'Time (ms)')
    #ax.set_xscale('log')
    plt.title(name)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Luminosity (ergs/s)')
    plt.show()

    percentage_events = events_count / len(Earth_Lum) * 100
    print(percentage_events, '% events above threshold')
    return percentage_events / 100

#%%

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

def plotlog(energy, E_v, X_ref, X_E, S_d, S_ref, L_ref, Earth_Lum, name):
    """
    

    Parameters
    ----------
    energy : array
        reference energy used to plot threshold curve
    E_v : array
        neutrino energy for the event in question
    X_ref : float
        reference cross section used to calculate threshold curve
    X_E : array
        cross section of neutrinos at the energies of the E_v array
    S_d : float
        surface area of event
    S_ref : float
        surface area of detector
    L_ref : float
        reference luminosity used to calculate threshold curve
    Earth_Lum : array
        neutrino luminosities of event in question at Earth 
        ie scaled depending on geometry
    name : string
        list element for title of graphs
    t : array
        times of events
    xlim : tuple
        domain of graph
    ylim : tuple
        range of graph 

    Returns
    -------
    events_count : float
        number of events above threshold

    """
    #find the threshold luminosities for our dataset
    L_t = lum_thresh(X_ref, X_E, S_d, S_ref, L_ref)
    #threshold is min L you see a few neutrinos
    
    logE = np.log(energy)
    logL_t = np.log(L_t)

    yerr = 10e45


    #calc number of events above threshold
    events_count = thresholdlog(Earth_Lum, L_t)


    #plotting the graph with all data points and threshold curve
    fig, ax = plt.subplots(figsize=(8, 8))
    plot = plt.scatter(E_v, Earth_Lum, color = 'magenta')
    plt.errorbar(energy, L_t, label = 'threshold', fmt = 'x', color = 'turquoise')
    plt.plot(energy, L_t, color = 'turquoise')
    plt.legend()
    ax.set_yscale('log')
    plt.semilogx()
    plt.title(name)
    plt.xlabel('Energy (MeV)')
    plt.ylabel('Luminosity (ergs/s)')
    plt.show()
    
    percentage_events = events_count / len(Earth_Lum) * 100
    print(percentage_events, '% events above threshold')
    return percentage_events / 100


def plot2(energy, E_v_i, E_v_c, X_ref, X_E, S_d, S_ref, L_ref, Earth_Lum_i, Earth_Lum_c, name, t_i, t_c, xlim, ylim, p0):
    """
    

    Parameters
    ----------
    energy : array
        reference energy used to plot threshold curve
    E_v_i : array
        neutrino energy for the irrotational event in question
    E_v_c : array
        neutrino energy for the corotational event in question
    X_ref :  float
        reference cross section used to calculate threshold curve
    X_E : array
        cross section of neutrinos at the energies of the E_v_i or E_v_c array
    S_d : float
        surface area of event
    S_ref :float
        surface area of detector
    L_ref : float
        reference luminosity used to calculate threshold curve
    Earth_Lum_i : array
        neutrino luminosities of irrotational event in question at Earth 
    Earth_Lum_c : array
        neutrino luminosities of corotational event in question at Earth 
    name : tuple
        list element for title of graphs
    t_i : array
        times of irrotational events
    t_c : array
        times of corotational events
    xlim : tuple
        domain of graph
    ylim : tuple
        range of graph 
        
    Returns
    -------
    TYPE: tuple
        events count for both types of event, [irrotational, corotational]

    """
    #find the threshold luminosities for our dataset
    L_t = lum_thresh(X_ref, X_E, S_d, S_ref, L_ref)
    #threshold is min L you see a few neutrinos

    #calc the piecewise fit for threshold
    p0 = p0
    param , cov = curve_fit(piecewise_linear, energy, L_t, p0)
    yerr = 10e50
    
    #calc number of events above threshold
    events_count_i = threshold(Earth_Lum_i, E_v_i, param)
    events_count_c = threshold(Earth_Lum_c, E_v_c, param)
    
    #plotting the graph with all data points and threshold curve
    fig, ax = plt.subplots(figsize=(14, 8))
    plot_i = plt.scatter(E_v_i, Earth_Lum_i, c = t_i, cmap=cm.get_cmap('spring'))
    plot_c = plt.scatter(E_v_c, Earth_Lum_c, c = t_c, cmap=cm.get_cmap('winter'))
    plt.errorbar(energy, L_t, label = 'threshold', yerr = yerr, fmt = 'x', color = 'turquoise')
    plt.plot(energy, piecewise_linear(energy, *param), color = 'turquoise')
    plt.legend()
    cb = plt.colorbar(plot_i, ax=ax, label = 'Time Irrotational (ms)')
    cb2 = plt.colorbar(plot_c, ax=ax, label = 'Time Corotational (ms)')
    ax.set_yscale('log')
    plt.title(name)
    #plt.xlim(xlim)
    #plt.ylim(ylim)
    plt.xlabel('Energy (MeV)')
    plt.ylabel( 'Luminosity (ergs/s)')
    plt.show()
    
    #print how many events are above threshold for irro and coro
    percentage_events_i = events_count_i / len(Earth_Lum_i) * 100
    print(percentage_events_i, '% events above threshold for irrotational')
    percentage_events_c = events_count_c / len(Earth_Lum_c) * 100
    print(percentage_events_c, '% events above threshold for corotational')
    return np.array([events_count_i, events_count_c])

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
        
        #plt.plot(pdf_x, pdf_y, label = 'Gaussian Best Fit', color = 'magenta')
        #plot histogram of energies
        plt.hist(energy_list, density=True, bins = 50, alpha = 0.25, color = 'turquoise', label = 'Raw Data')  
        plt.ylabel('Events')
        plt.title(name +' Histogram')
        plt.xlabel('Energies (MeV)')
        plt.legend()
        plt.show()
        
        return energy_list
    

#%% determining references for scale


def events(X_ref, L_ref, t_ref, no_events_ref, X_E, L_Earth_E, E, percent, name):
    
    cf = conversion_factor(X_ref, L_ref, t_ref, no_events_ref)
    actual_events = np.array(rate_events(X_E, L_Earth_E)) * cf * percent

    #prop_list = proportion_per_bin(actual_events)
    
    plt.title('Histogram 1')
    hist, bins, idk = plt.hist(E, bins = len(E), weights = actual_events)
    plt.show()
    
    
    histactual =  []
    for i in range(0, len(E)):
        histactual.append(int(hist[i]) * actual_events[i])
    
    s2 = np.random.poisson(lam=(actual_events), size = len(E))
    
    plt.title('Histogram for Poisson Fluctuations for '  
              + name)
    hist2, bins, idk = plt.hist(E, weights = (s2), bins = 10, color = 'turquoise', alpha = 0.25, label = 'Events')
    plt.xlabel('Energy (MeV)')
    plt.legend()
    plt.ylabel('Events')
    plt.show()


    
def proportion_per_bin(onedlist):
    total_events = 0
    for i in range(0, len(onedlist)):
        total_events =+ onedlist[i]

    prop_list = []
    for i in range(0, len(onedlist)):
        prop_list.append(int((onedlist[i]/total_events) * 100))
        
    prop_list = np.array(prop_list)
    #print(prop_list)
    return prop_list


#%%
def distgraph(energy, E_v, X_ref, X_E, S_d, S_ref, L_ref, Earth_Lum, p0):
    """
    

    Parameters
    ----------
    energy : array
        reference energy used to plot threshold curve
    E_v : array
        neutrino energy for the event in question
    X_ref : float
        reference cross section used to calculate threshold curve
    X_E : array
        cross section of neutrinos at the energies of the E_v array
    S_d : float
        surface area of event
    S_ref : float
        surface area of detector
    L_ref : float
        reference luminosity used to calculate threshold curve
    Earth_Lum : array
        neutrino luminosities of event in question at Earth 
        ie scaled depending on geometry


    Returns
    -------
    events_count : float
        number of events above threshold

    """
    #find the threshold luminosities for our dataset
    L_t = lum_thresh(X_ref, X_E, S_d, S_ref, L_ref)
    #threshold is min L you see a few neutrinos
    
    
    p0 = p0
    param , cov = curve_fit(piecewise_linear, energy, L_t, p0)
    yerr = 10e50
    
    #calc number of events above threshold
    events_count = threshold(Earth_Lum, E_v, param)
    
    percentage_events = events_count / len(Earth_Lum) * 100
    print(percentage_events, '% events above threshold')
    return percentage_events / 100






"""
Created on Thu Jul  7 13:38:23 2022

@author: trinitystenhouse

Here is an example of what the code can be used to do. This was for a research 
project in the 2nd year of my 4 year integrated Masters degree, called a UROP.
I read in data from various merger events and then plotted graphs and did simulations
to get energy spectra at earth, max distances we could see neutrinos until and
investigate angular dependence. You, of course, can do whatever you like with it.
Please bear in mind I was never taught how to code and thhis was all improvised, 
so if it isn't really efficient I'm sorry, and feel free to send me recommendations.
"""

import random
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from Geometries2 import *
from luminosities2 import *


#%% BNS merger data

#irrotational BNS Dataset
t_BNS_i, L_ve_BNS_i, E_ve_BNS_i, X_ve_BNS_i, L_vae_BNS_i, E_vae_BNS_i, X_vae_BNS_i, \
    L_vx_BNS_i, E_vx_BNS_i, X_vx_BNS_i = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BNS_data_irro.txt', skiprows = 1, unpack = True)
#corotational BNS Dataset
t_BNS_c, L_ve_BNS_c, E_ve_BNS_c, X_ve_BNS_c, L_vae_BNS_c, E_vae_BNS_c, X_vae_BNS_c, \
    L_vx_BNS_c, E_vx_BNS_c, X_vx_BNS_c = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BNS_data_coro2.txt', skiprows = 1, unpack = True)
#irrotational BNS Reference Dataset
energy, Lve, Xve, Lvae, Xvae, Lvx, Xvx = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/energybins.txt', skiprows = 1, unpack = True)
#irrotational BNS Reference Dataset
energy2, Lve2, Xve2, Lvae2, Xvae2, Lvx2, Xvx2 = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/energybinscoro.txt', skiprows = 1, unpack = True)

#input distance to event kpc, angle from plane of events
Geo_BNS = Torus(10)  #ref geometry
S_BNS = Geo_BNS.Area()

L_Earth_ve_BNS_i = Geo_BNS.Earth_Lum(L_ve_BNS_i * 1000)
L_Earth_vae_BNS_i = Geo_BNS.Earth_Lum(L_vae_BNS_i * 1000)
L_Earth_vx_BNS_i = Geo_BNS.Earth_Lum(L_vx_BNS_i * 1000)
L_Earth_ve_BNS_c = Geo_BNS.Earth_Lum(L_ve_BNS_c * 1000)
L_Earth_vae_BNS_c = Geo_BNS.Earth_Lum(L_vae_BNS_c * 1000)
L_Earth_vx_BNS_c = Geo_BNS.Earth_Lum(L_vx_BNS_c * 1000)
#luminosity in ergs/s

#converting cross section to correct powers
X_ve_i = X_ve_BNS_i * 10 ** -38
X_vae_i = X_vae_BNS_i * 10 ** -38
X_vx_i = X_vx_BNS_i * 10 ** -38 #all in cm^2
X_ve_c = X_ve_BNS_c * 10 ** -38
X_vae_c = X_vae_BNS_c * 10 ** -38
X_vx_c = X_vx_BNS_c * 10 ** -38 #all in cm^2

Xve = Xve * 10 ** -38
Xvae = Xvae * 10 ** -38
Xvx = Xvx * 10 ** -38
 
Xve2 = Xve2 * 10 ** -38
Xvae2 = Xvae2 * 10 ** -38
Xvx2 = Xvx2 * 10 ** -38

namelisti = ['Irrotational \n Electron Neutrino', 'Irrotational \n Electron AntiNeutrino', 'Irrotational \n Non Electron Neutrino']
namelistc = ['Corotational \n Electron Neutrino', 'Corotational \n Electron AntiNeutrino', 'Corotational \n Non Electron Neutrino']

#%% BBH merger data

E_GeV_BBH, L_BBH_1_opt, L_BBH_1_un, X, L_BBH_3_opt, L_BBH_3_un = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BBH_data.txt', skiprows = 1, unpack = True)
#luminosity in ergs/s

E_MeV_BBH = E_GeV_BBH * 1000

#%% Reference Supernova data
#using DUNE electron-capture supernova neutrino luminosities and data
            #as reference values
Sphere_ref = Spherical(100) #using supernova at 100kpc
S_ref = Sphere_ref.Area()
L_ref = 0.5 * 10 ** 52
Earth_Lum_ref = max(Sphere_ref.Earth_Lum([L_ref]))
E_ref = 9.5
#10 interactions in 10s

X_ve_ref = 10 ** -5 * 10 ** -38
X_vae_ref = 13 ** -6 * 10 ** -38
X_vx_ref = 10 ** -6 * 10 ** -38
#%% BNS Reference data
refdataBNS_ve = [energy, np.array(X_ve_ref), 100, L_ref]
refdataBNS_vae = [energy, np.array(X_vae_ref), 8, L_ref]
refdataBNS_vx = [energy, np.array(X_vx_ref), 40, L_ref]

eventsdata_ve_i = [E_ve_BNS_i, Xve, S_BNS, L_Earth_ve_BNS_i]
eventsdata_ve_c = [E_ve_BNS_c, Xve2, S_BNS, L_Earth_ve_BNS_c]
eventsdata_vae_i = [E_vae_BNS_i, Xvae, S_BNS, L_Earth_vae_BNS_i]
eventsdata_vae_c =[E_vae_BNS_c, Xvae2, S_BNS, L_Earth_vae_BNS_c]
eventsdata_vx_i = [E_vx_BNS_i, Xvx, S_BNS, L_Earth_vx_BNS_i]
eventsdata_vx_c = [E_vx_BNS_c, Xvx2, S_BNS, L_Earth_vx_BNS_c]
full_events_ve = [E_ve_BNS_i, Xve, S_BNS, L_Earth_ve_BNS_i, E_ve_BNS_c, Xve2, S_BNS, L_Earth_ve_BNS_c]
full_events_vae = [E_vae_BNS_i, Xvae, S_BNS, L_Earth_vae_BNS_i, E_vae_BNS_c, Xvae2, S_BNS, L_Earth_vae_BNS_c]
full_events_vx = [E_vx_BNS_i, Xvx, S_BNS, L_Earth_vx_BNS_i, E_vx_BNS_c, Xvx2, S_BNS, L_Earth_vx_BNS_c]

p0 = [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]
#%% Electron Neutrino
events_count_ve_i = full(refdataBNS_ve, eventsdata_ve_i,['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise')
events_count_ve_c = full(refdataBNS_ve, eventsdata_ve_c, ['Corotational BNS Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 4.5e52], p0], 'piecewise')

array_ve = full(refdataBNS_ve, full_events_ve, ['BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 4.5e52], p0, t_BNS_c,], 'multi')

#Poisson Fluctuations
events(refdataBNS_ve, eventsdata_ve_i, np.array(Earth_Lum_ref), X_ve_i, 10, 10, events_count_ve_i, namelisti[0])
events(refdataBNS_ve, eventsdata_ve_c, np.array(Earth_Lum_ref), X_ve_c, 10, 10, events_count_ve_c, namelistc[0])

#Monte Carlo Sim
#note first number should be small if high L since otherwise uses too much computational power
#energy_hist(1e-7, X_ve_ref, X_ve_i, L_ref, L_Earth_ve_BNS_i, 10, 10, events_count_ve_i, E_ve_BNS_i, 0.3, namelisti[0])
#energy_hist(1e-7, X_ve_ref, X_ve_c, L_ref, L_Earth_ve_BNS_c, 10, 10, events_count_ve_c, E_ve_BNS_c, 0.3, namelistc[0])

#%% Electron AntiNeutrino
events_count_vae_i = full(refdataBNS_vae, eventsdata_vae_i,['Irrotational BNS Electron AntiNeutrino', t_BNS_i, [0.5, 30], [-0.1e52, 8e52], p0], 'piecewise')
events_count_vae_c = full(refdataBNS_vae, eventsdata_vae_c, ['Corotational BNS Electron AntiNeutrino', t_BNS_c, [0.5, 30], [-0.1e52, 8e52], p0], 'piecewise')

array_vae = full(refdataBNS_vae, full_events_vae, ['BNS Electron AntiNeutrino', t_BNS_i, [0.5, 30], [-0.1e52, 8e52], p0, t_BNS_c,], 'multi')

#Poisson Fluctuations
events(refdataBNS_vae, eventsdata_vae_i, np.array(Earth_Lum_ref), X_vae_i, 10, 10, events_count_vae_i, namelisti[1])
events(refdataBNS_vae, eventsdata_vae_c, np.array(Earth_Lum_ref), X_vae_c, 10, 10, events_count_vae_c, namelistc[1])

#Monte Carlo Sim
#energy_hist(1e-7, X_vae_ref, X_vae_i, L_ref, L_Earth_vae_BNS_i, 10, 10, events_count_vae_i, E_vae_BNS_i, 0.27, namelisti[1])
#energy_hist(1e-7, X_vae_ref, X_vae_c, L_ref, L_Earth_vae_BNS_c, 10, 10, events_count_vae_c, E_vae_BNS_c, 0.3, namelistc[1])


#%% Non-Electron Neutrino 
events_count_vx_i = full(refdataBNS_vx, eventsdata_vx_i,['Irrotational BNS NonElectron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 4.5e52], p0], 'piecewise')
events_count_vx_c = full(refdataBNS_vx, eventsdata_vx_c, ['Corotational BNS NonElectron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise')

array_vx = full(refdataBNS_vx, full_events_vx, ['BNS NonElectron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 4.5e52], p0, t_BNS_c,], 'multi')

#Poisson Fluctuations
events(refdataBNS_vx, eventsdata_vx_i, np.array(Earth_Lum_ref), X_vx_i, 10, 10, events_count_vx_i, namelisti[2])
events(refdataBNS_vx, eventsdata_vx_c, np.array(Earth_Lum_ref), X_vx_c, 10, 10, events_count_vx_c, namelistc[2])

#Monte Carlo Sim
#energy_hist(1e-7, X_vx_ref, X_vx_i, L_ref, L_Earth_vx_BNS_i, 10, 10, events_count_vx_i, E_vx_BNS_i, 0.2, namelisti[2])
#energy_hist(1e-7, X_vx_ref, X_vx_c, L_ref, L_Earth_vx_BNS_c, 10, 10, events_count_vx_c, E_vx_BNS_c, 0.3, namelistc[2])

#%% BBH Data

Geo_BBH = Torus(1, 0)
S_BBH = Geo_BBH.Area()
L_Earth_BBH_3_opt = Geo_BBH.Earth_Lum(L_BBH_3_opt)
L_Earth_BBH_1_opt = Geo_BBH.Earth_Lum(L_BBH_1_opt)
L_Earth_BBH_3_un = Geo_BBH.Earth_Lum(L_BBH_3_un)
L_Earth_BBH_1_un = Geo_BBH.Earth_Lum(L_BBH_1_un)


refdataBBH = [E_MeV_BBH, X_ve_ref, 814.8, L_ref]
eventdata_BBH1opt = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_1_opt]
eventdata_BBH1un = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_1_un]
eventdata_BBH3opt = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_3_opt]
eventdata_BBH3un = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_3_un]
#%% Black Hole Merger 3x10^7 Msun

count_BBH3_opt = full(refdataBBH, eventdata_BBH3opt, ['$3x10^{7}$ M$_{sun}$ BBH (optimistic)'], 'log')
count_BBH3_un = full(refdataBBH, eventdata_BBH3un, ['$3x10^7$ M$_{sun}$ BBH (unoptimistic)'], 'log')

#Poisson Fluctuations
events(refdataBBH, eventdata_BBH3opt, np.array(Earth_Lum_ref), X, 10, 10, count_BBH3_opt, '$3x10^{7}$ M$_{sun}$ \n BBH (optimistic)')
events(refdataBBH, eventdata_BBH3un, np.array(Earth_Lum_ref), X, 10, 10, count_BBH3_un, '$3x10^7$ M$_{sun}$ \n BBH (unoptimistic)')

#%%Black Hole Merger 1x10^7 Msun

count_BBH1_opt = full(refdataBBH, eventdata_BBH1opt, ['$1x10^7$ M$_{sun}$ BBH (optimistic)'], 'log')
count_BBH1_un = full(refdataBBH, eventdata_BBH1un, ['$1x10^7$ M$_{sun}$ BBH (unoptimistic)'], 'log')

#Poisson Fluctuations 
events(refdataBBH, eventdata_BBH1opt, np.array(Earth_Lum_ref), X, 10, 10, count_BBH1_opt, '$1x10^7$ M$_{sun}$ \n BBH (optimistic)')
events(refdataBBH, eventdata_BBH1un, np.array(Earth_Lum_ref), X, 10, 10, count_BBH1_un, '$1x10^7$ M$_{sun}$ \n BBH (unoptimistic)')

#%%For comparison at different angles
refdataBBH = [E_MeV_BBH, X_ve_ref, 814.8, L_ref]

Geo_BBH_pi2 = Torus(10, (89/360 * 2 * np.pi))
S_BBH_pi2 = Geo_BBH_pi2.Area()
L_Earth_BBH_3_opt_pi2 = Geo_BBH_pi2.Earth_Lum(L_BBH_3_opt)
L_Earth_BBH_1_opt_pi2 = Geo_BBH_pi2.Earth_Lum(L_BBH_1_opt)
eventdata_BBH1opt_angle = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_1_opt, E_MeV_BBH, X, S_BBH_pi2, L_Earth_BBH_1_opt_pi2]
eventdata_BBH3opt_angle = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_3_opt, E_MeV_BBH, X, S_BBH_pi2, L_Earth_BBH_3_opt_pi2]
L_Earth_BBH_3_un_pi2 = Geo_BBH_pi2.Earth_Lum(L_BBH_3_un)
L_Earth_BBH_1_un_pi2 = Geo_BBH_pi2.Earth_Lum(L_BBH_1_un)
eventdata_BBH1un_angle = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_1_un, E_MeV_BBH, X, S_BBH_pi2, L_Earth_BBH_1_un_pi2]
eventdata_BBH3un_angle = [E_MeV_BBH, X, S_BBH, L_Earth_BBH_3_un, E_MeV_BBH, X, S_BBH_pi2, L_Earth_BBH_3_un_pi2]

count_BBH1_opt_angles = full(refdataBBH, eventdata_BBH1opt_angle, ['$1x10^7$ M$_{sun}$ BBH (optimistic)'], 'multilog')
count_BBH3_opt_angles = full(refdataBBH, eventdata_BBH3opt_angle, ['$3x10^7$ M$_{sun}$ BBH (optimistic)'], 'multilog')
count_BBH1_un_angles = full(refdataBBH, eventdata_BBH1un_angle, ['$1x10^7$ M$_{sun}$ BBH (unoptimistic)'], 'multilog')
count_BBH3_un_angles = full(refdataBBH, eventdata_BBH3un_angle, ['$3x10^7$ M$_{sun}$ BBH (unoptimistic)'], 'multilog')


#%% generate random dataset

RP = random.randint(50,60)
RN = random.randint(1,100) * 10 ** RP
print(RN)
Space_Lum_1 = ( L_ve_BNS_i / 10 ** 52 ) * RN #luminosity of an event at the event
RD = random.randint(1,20)
Sphere_Geo = Spherical(RD) #input distance to event
S_d = Sphere_Geo.Area() #ref geometry
Earth_Lum_1 = Sphere_Geo.Earth_Lum(Space_Lum_1)

X_E_1 = []
for i in range(0,len(Xve)):
    n = random.randint(1,20) * 10 ** -40
    X_E_1.append(n)
    

E_1 = []
for i in range(0,len(E_ve_BNS_i)):
    n = random.randint(1,300) /10
    E_1.append(n)
    
if RD < 10:
    print('Random distance to event ', RD, 'kpc, less than ref')
else:
    print('Random distance to event ', RD, 'kpc, greater than ref')
    
if np.median(X_E_1) < np.median(X_ve_i):
    print('Random Cross Sections ', X_E_1, 'in cm^2, less than')
else:
    print('Random Cross Sections ', X_E_1, 'in cm^2, more than')
    
if np.median(Space_Lum_1) < np.median(L_ve_BNS_i):   
    print('Random Luminosities ', Space_Lum_1, 'ergs/s, less than')
else:
    print('Random Luminosities ', Space_Lum_1, 'ergs/s, more than')
    
if np.median(E_1) < np.median(E_ve_BNS_i):   
    print('Random Luminosities ', E_1, 'ergs/s, less than')
else:
    print('Random Luminosities ', E_1, 'ergs/s, more than')

randevents = [E_1, X_E_1, S_d, Earth_Lum_1]
#%% Random Dataset Plots
events_count_rand = full(refdataBNS_ve, randevents, ['Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 1e53], p0], 'piecewise')
events(refdataBNS_ve, randevents, np.array(Earth_Lum_ref), X_ve_i, 10, 10, events_count_rand, 'Randomised Electron Neutrino')

events_count_rand2 = full(refdataBNS_vae, randevents, ['Electron AntiNeutrino', t_BNS_c, [0.5, 30], [-0.1e52, 1e53], p0], 'piecewise')
events(refdataBNS_vae, randevents, np.array(Earth_Lum_ref), X_vae_i, 10, 10, events_count_rand2, 'Randomised Electron AntiNeutrino')


#%% looking at the maximum distance the merger can be to still see neutrinos - BNS

distancelist = np.arange(0.0,100.0,10.0)

list_ve_i = []
list_ve_c = []
list_vae_i = []
list_vae_c = []
list_vx_i = []
list_vx_c = []

for i in range(distancelist.shape[0]):
    geo = Torus(distancelist[i], 0) #can change the angle
    area = geo.Area()
    list_ve_i.append(distgraph(refdataBNS_ve, eventsdata_ve_i, area, geo.Earth_Lum(L_ve_BNS_i * 1000), 'Irrotational ve', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_ve_c.append(distgraph(refdataBNS_ve, eventsdata_ve_c, area, geo.Earth_Lum(L_ve_BNS_c * 1000), 'Corotational ve', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))                                                                                                                                                                                        
    list_vae_i.append(distgraph(refdataBNS_vae, eventsdata_vae_i, area, geo.Earth_Lum(L_vae_BNS_i * 1000),'Irrotational vae', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_vae_c.append(distgraph(refdataBNS_vae, eventsdata_vae_c, area, geo.Earth_Lum(L_vae_BNS_c * 1000), 'Cotational vae', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_vx_i.append(distgraph(refdataBNS_vx, eventsdata_vx_i, area, geo.Earth_Lum(L_vx_BNS_i * 1000), 'Irrotational vx', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_vx_c.append(distgraph(refdataBNS_vx,  eventsdata_vx_c, area, geo.Earth_Lum(L_vx_BNS_c * 1000),'Cotational vx', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))

plt.plot(distancelist, list_ve_i, marker = 'x', label = 'Irrotational ve', color = 'magenta')
plt.plot(distancelist, list_ve_c, marker = 'x', label = 'Corotational ve', color = 'violet')
plt.plot(distancelist, list_vae_i, marker = 'x', label = 'Irrotational vae', color = 'aqua')
plt.plot(distancelist, list_vae_c, marker = 'x', label = 'Corotational vae', color = 'turquoise')
plt.plot(distancelist, list_vx_i, marker = 'x', label = 'Irrotational vx', color = 'darkviolet')
plt.plot(distancelist, list_vx_c, marker = 'x', label = 'Corotational vx', color = 'mediumpurple')
plt.legend()
plt.xlabel('Distance of Merger (kpc)')
plt.ylabel('Proportion of Events Above Threshold')
plt.show()

#%% angular dependence of BNS merger
angleslist = np.arange(0.0,np.pi/2,np.pi/40.0)

list_ve_i = []
list_ve_c = []
list_vae_i = []
list_vae_c = []
list_vx_i = []
list_vx_c = []

for i in range(angleslist.shape[0]):
    geo = Torus(50, angleslist[i]) #can change the distance
    area = geo.Area()
    list_ve_i.append(distgraph(refdataBNS_ve, eventsdata_ve_i, area, geo.Earth_Lum(L_ve_BNS_i * 1000), 'Irrotational ve', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_ve_c.append(distgraph(refdataBNS_ve, eventsdata_ve_c, area, geo.Earth_Lum(L_ve_BNS_c * 1000), 'Corotational ve', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))                                                                                                                                                                                        
    list_vae_i.append(distgraph(refdataBNS_vae, eventsdata_vae_i, area, geo.Earth_Lum(L_vae_BNS_i * 1000),'Irrotational vae', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_vae_c.append(distgraph(refdataBNS_vae, eventsdata_vae_c, area, geo.Earth_Lum(L_vae_BNS_c * 1000), 'Cotational vae', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_vx_i.append(distgraph(refdataBNS_vx, eventsdata_vx_i, area, geo.Earth_Lum(L_vx_BNS_i * 1000), 'Irrotational vx', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))
    list_vx_c.append(distgraph(refdataBNS_vx,  eventsdata_vx_c, area, geo.Earth_Lum(L_vx_BNS_c * 1000),'Cotational vx', ['Irrotational BNS Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], p0], 'piecewise'))

plt.plot(angleslist, list_ve_i, marker = 'x', label = 'Irrotational ve', color = 'magenta')
plt.plot(angleslist, list_ve_c, marker = 'x', label = 'Corotational ve', color = 'violet')
plt.plot(angleslist, list_vae_i, marker = 'x', label = 'Irrotational vae', color = 'aqua')
plt.plot(angleslist, list_vae_c, marker = 'x', label = 'Corotational vae', color = 'turquoise')
plt.plot(angleslist, list_vx_i, marker = 'x', label = 'Irrotational vx', color = 'darkviolet')
plt.plot(angleslist, list_vx_c, marker = 'x', label = 'Corotational vx', color = 'mediumpurple')
plt.legend()
plt.xlabel('Angles of Earth from Merger Plane (rad)')
plt.ylabel('Proportion of Events Above Threshold')
plt.show()

#%% distance we can see the merger until - BBH
distancelist = np.arange(0.0,10.0,1.0)
print(distancelist)
list_BBH3opt = []
list_BBH3un = []
list_BBH1opt = []
list_BBH1un = []


for i in range(distancelist.shape[0]):
    geo = Torus(distancelist[i], 0) #can change angle
    area = geo.Area()
    list_BBH3opt.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_3_opt * 1000), 'name', [0], 'log'))
    list_BBH3un.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_3_un * 1000), 'name', [0], 'log'))                                                                                                                                                                                          
    list_BBH1opt.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_1_opt * 1000), 'name', [0],'log'))
    list_BBH1un.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_1_un * 1000),'name',  [0], 'log'))

plt.plot(distancelist, list_BBH3opt, marker = 'x', label = 'Optimistic 3x10^7 mass BBH', color = 'magenta')
plt.plot(distancelist, list_BBH3un, marker = 'x', label = 'Unoptimistic 3x10^7 mass BBH', color = 'violet')
plt.plot(distancelist, list_BBH1opt, marker = 'x', label = 'Optimistic 1x10^7 mass BBH', color = 'aqua')
plt.plot(distancelist, list_BBH1un, marker = 'x', label = 'Unoptimistic 1x10^7 mass BBH', color = 'turquoise')
plt.legend()
plt.xlabel('Distance of Merger (kpc)')
plt.ylabel('Proportion of Events Above Threshold')
plt.show()

#%% angular dependence of BBH
angleslist = np.arange(0.0,np.pi/2,np.pi/40.0)

list_BBH3opt = []
list_BBH3un = []
list_BBH1opt = []
list_BBH1un = []

for i in range(angleslist.shape[0]):
    geo = Torus(1, angleslist[i]) #can change distance
    area = geo.Area()
    list_BBH3opt.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_3_opt * 1000), 'name', [0], 'log'))
    list_BBH3un.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_3_un * 1000), 'name', [0], 'log'))                                                                                                                                                                                          
    list_BBH1opt.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_1_opt * 1000), 'name', [0],'log'))
    list_BBH1un.append(distgraph(refdataBBH, eventdata_BBH3opt, area, geo.Earth_Lum(L_BBH_1_un * 1000),'name',  [0], 'log'))


plt.plot(angleslist, list_BBH3opt, marker = 'x', label = 'Optimistic 3x10^7 mass BBH', color = 'magenta')
plt.plot(angleslist, list_BBH3un, marker = 'x', label = 'Unoptimistic 3x10^7 mass BBH', color = 'violet')
plt.plot(angleslist, list_BBH1opt, marker = 'x', label = 'Optimistic 1x10^7 mass BBH', color = 'aqua')
plt.plot(angleslist, list_BBH1un, marker = 'x', label = 'Unoptimistic 1x10^7 mass BBH', color = 'turquoise')
plt.legend()
plt.xlabel('Angles of Earth from Merger Plane (rad)')
plt.ylabel('Proportion of Events Above Threshold')
plt.show()

#%% checking conversion factor works
cf = conversion_factor(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10)
rate = np.array(X_ve_ref) * Earth_Lum_ref
g = cf * rate * 10
print(g)

plt.bar(9.5, g, color='turquoise', width = 0.1)
plt.show()

#%% checking torus class has correct dependencies
distances1 = np.arange(0.0,500.0,50.0)
angles1 = np.arange(0.0,np.pi/2,np.pi/20.0)
the_array1 = np.zeros((angles1.shape[0], distances1.shape[0]))

for i in range(distances1.shape[0]):
    for j in range(angles1.shape[0]):
   
        TorusA = Torus(distances1[i], angles1[j])
        the_array1[j][i] = TorusA.Earth_Lum(([1]))[0] 
        print(distances1[i], angles1[j], the_array1[i][j])
      
plt.pcolormesh(distances1, angles1, the_array1)
plt.colorbar()
plt.title('Torus Colour Plot')
plt.xlabel('distance (kpc)')
plt.ylabel('angle (radians)')
plt.show()
#%% checking spherical class has correct dependencies
distances2 = np.arange(0.0,1000.0,100.0)
angles2 = np.arange(0.0,np.pi/2,np.pi/200.0)
the_array2 = np.zeros((angles2.shape[0], distances2.shape[0]))

for i in range(distances2.shape[0]):
    for j in range(angles2.shape[0]):
   
        SphereA = Spherical(distances2[i])
        the_array2[j][i] = SphereA.Earth_Lum(([1]))[0] 
        print(distances2[i], angles2[j], the_array2[j][i])

#try without function just use the numbers

#SphereA.Earth_Lum(([1]))[0]
plt.pcolormesh(distances2, angles2, the_array2)
plt.colorbar()
plt.title('Sphere Colour Plot')
plt.xlabel('distance (kpc)')
plt.ylabel('angle (radians)')


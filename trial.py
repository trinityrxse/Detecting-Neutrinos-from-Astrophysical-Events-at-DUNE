
import random
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from Geometries2 import *
from luminosities2 import *


#%%
#using DUNE electron-capture supernova neutrino luminosities and data
            #as reference values

t_BNS_i, L_ve_BNS_i, E_ve_BNS_i, X_ve_BNS_i, L_vae_BNS_i, E_vae_BNS_i, X_vae_BNS_i, \
    L_vx_BNS_i, E_vx_BNS_i, X_vx_BNS_i = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BNS_data_irro.txt', skiprows = 1, unpack = True)

t_BNS_c, L_ve_BNS_c, E_ve_BNS_c, X_ve_BNS_c, L_vae_BNS_c, E_vae_BNS_c, X_vae_BNS_c, \
    L_vx_BNS_c, E_vx_BNS_c, X_vx_BNS_c = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BNS_data_coro2.txt', skiprows = 1, unpack = True)

E_GeV_BBH, L_BBH_1_opt, L_BBH_1_un, X, L_BBH_3_opt, L_BBH_3_un = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BBH_data.txt', skiprows = 1, unpack = True)

energy, Lve, Xve, Lvae, Xvae, Lvx, Xvx = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/energybins.txt', skiprows = 1, unpack = True)
energy2, Lve2, Xve2, Lvae2, Xvae2, Lvx2, Xvx2 = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/energybinscoro.txt', skiprows = 1, unpack = True)

#luminosity in ergs/s

Xve = Xve * 10 ** -38
Xvae = Xvae * 10 ** -38
Xvx = Xvx * 10 ** -38
 
Xve2 = Xve2 * 10 ** -38
Xvae2 = Xvae2 * 10 ** -38
Xvx2 = Xvx2 * 10 ** -38

#input distance to event kpc, angle from plane of events
Geo_BNS = Torus(10)  #ref geometry
S_BNS = Geo_BNS.Area()

L_Earth_ve_BNS_i = Geo_BNS.Earth_Lum(L_ve_BNS_i * 1000)
L_Earth_vae_BNS_i = Geo_BNS.Earth_Lum(L_vae_BNS_i * 1000)
L_Earth_vx_BNS_i = Geo_BNS.Earth_Lum(L_vx_BNS_i * 1000)
L_Earth_ve_BNS_c = Geo_BNS.Earth_Lum(L_ve_BNS_c * 1000)
L_Earth_vae_BNS_c = Geo_BNS.Earth_Lum(L_vae_BNS_c * 1000)
L_Earth_vx_BNS_c = Geo_BNS.Earth_Lum(L_vx_BNS_c * 1000)


X_ve_i = X_ve_BNS_i * 10 ** -38
X_vae_i = X_vae_BNS_i * 10 ** -38
X_vx_i = X_vx_BNS_i * 10 ** -38 #all in cm^2
X_ve_c = X_ve_BNS_c * 10 ** -38
X_vae_c = X_vae_BNS_c * 10 ** -38
X_vx_c = X_vx_BNS_c * 10 ** -38 #all in cm^2

namelisti = ['Irrotational \n Electron Neutrino', 'Irrotational \n Electron AntiNeutrino', 'Irrotational \n Non Electron Neutrino']

namelistc = ['Corotational \n Electron Neutrino', 'Corotational \n Electron AntiNeutrino', 'Corotational \n Non Electron Neutrino']

Geo_BBH = Torus(0.5, 1.4)
S_BBH = Geo_BBH.Area()
L_Earth_BBH_3_opt = Geo_BBH.Earth_Lum(L_BBH_3_opt)
L_Earth_BBH_1_opt = Geo_BBH.Earth_Lum(L_BBH_1_opt)
L_Earth_BBH_3_un = Geo_BBH.Earth_Lum(L_BBH_3_un)
L_Earth_BBH_1_un = Geo_BBH.Earth_Lum(L_BBH_1_un)

E_MeV_BBH = E_GeV_BBH * 1000

Sphere_ref = Spherical(100) #using supernova at 100kpc
S_ref = Sphere_ref.Area()
L_ref = 0.5 * 10 ** 52
Earth_Lum_ref = max(Sphere_ref.Earth_Lum([L_ref]))
E_ref = 9.5
#10 interactions in 10s

X_ve_ref = 10 ** -5 * 10 ** -38
X_vae_ref = 13 ** -6 * 10 ** -38
X_vx_ref = 10 ** -6 * 10 ** -38

#1.9445e-17
#%% Electron Neutrino
events_count_ve_i = plot(energy, E_ve_BNS_i, X_ve_ref, Xve, S_BNS, 814.8, L_ref, L_Earth_ve_BNS_i,'Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])
events_count_ve_c = plot(energy2, E_ve_BNS_c, X_ve_ref, Xve2, S_BNS, 814.8, L_ref, L_Earth_ve_BNS_c,'Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 4.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

array_ve = plot2(energy, E_ve_BNS_i, E_ve_BNS_c, X_ve_ref, Xve, S_BNS, 814.8, L_ref, L_Earth_ve_BNS_i, L_Earth_ve_BNS_c, 'Electron Neutrino', t_BNS_i, t_BNS_c, [0.5, 30], [-0.1e52, 4.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

#Poisson Fluctuations
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X_ve_i, L_Earth_ve_BNS_i, E_ve_BNS_i, events_count_ve_i, namelisti[0])
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X_ve_c, L_Earth_ve_BNS_c, E_ve_BNS_c, events_count_ve_c, namelistc[0])

#Monte Carlo Sim
#energy_hist(1e-7, X_ve_ref, X_ve_i, L_ref, L_Earth_ve_BNS_i, 10, 10, events_count_ve_i, E_ve_BNS_i, 0.3, namelisti[0])
#energy_hist(1e-7, X_ve_ref, X_ve_c, L_ref, L_Earth_ve_BNS_c, 10, 10, events_count_ve_c, E_ve_BNS_c, 0.3, namelistc[0])

#%% Electron AntiNeutrino
events_count_vae_i = plot(energy, E_vae_BNS_i, X_vae_ref, Xvae, S_BNS, 814.8, L_ref, L_Earth_vae_BNS_i, 'Electron AntiNeutrino', t_BNS_i, [0.5, 30], [-0.1e52, 8.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])
events_count_vae_c = plot(energy2, E_vae_BNS_c, X_vae_ref, Xvae2, S_BNS, 814.8, L_ref, L_Earth_vae_BNS_c, 'Electron AntiNeutrino', t_BNS_c, [0.5, 30], [-0.1e52, 8.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

array_vae = plot2(energy, E_vae_BNS_i, E_vae_BNS_c, X_vae_ref, Xvae, S_BNS, 814.8, L_ref, L_Earth_vae_BNS_i, L_Earth_vae_BNS_c, 'Electron AntiNeutrino', t_BNS_i, t_BNS_c, [0.5, 30], [-0.1e52, 8.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

#Poisson Fluctuations
events(np.array(X_vae_ref), np.array(Earth_Lum_ref), 10, 10, X_vae_i, L_Earth_vae_BNS_i, E_vae_BNS_i, events_count_vae_i, namelisti[1])
events(np.array(X_vae_ref), np.array(Earth_Lum_ref), 10, 10, X_vae_c, L_Earth_vae_BNS_c, E_vae_BNS_c, events_count_vae_c, namelistc[1])

#Monte Carlo Sim
#energy_hist(1e-7, X_vae_ref, X_vae_i, L_ref, L_Earth_vae_BNS_i, 10, 10, events_count_vae_i, E_vae_BNS_i, 0.27, namelisti[1])
#energy_hist(1e-7, X_vae_ref, X_vae_c, L_ref, L_Earth_vae_BNS_c, 10, 10, events_count_vae_c, E_vae_BNS_c, 0.3, namelistc[1])


#%% Non-Electron Neutrino 
events_count_vx_i = plot(energy, E_vx_BNS_i, X_vx_ref, Xvx, S_BNS, 814.8, L_ref, L_Earth_vx_BNS_i, 'Non-Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 5.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])
events_count_vx_c = plot(energy2, E_vx_BNS_c, X_vx_ref, Xvx2, S_BNS, 814.8, L_ref, L_Earth_vx_BNS_c, 'Non-Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 5.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

array_vx = plot2(energy, E_vx_BNS_i, E_vx_BNS_c, X_vx_ref, Xvx, S_BNS, 814.8, L_ref, L_Earth_vx_BNS_i, L_Earth_vx_BNS_c, 'Non-Electron Neutrino', t_BNS_i, t_BNS_c, [0.5, 30], [-0.1e52, 5.5e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

#Poisson Fluctuations
events(np.array(X_vx_ref), np.array(Earth_Lum_ref), 10, 10, X_vx_i, L_Earth_vx_BNS_i, E_vx_BNS_i, events_count_vx_i, namelisti[2])
events(np.array(X_vx_ref), np.array(Earth_Lum_ref), 10, 10, X_vx_c, L_Earth_vx_BNS_c, E_vx_BNS_c, events_count_vx_c, namelistc[2])

#Monte Carlo Sim
#energy_hist(1e-7, X_vx_ref, X_vx_i, L_ref, L_Earth_vx_BNS_i, 10, 10, events_count_vx_i, E_vx_BNS_i, 0.2, namelisti[2])
#energy_hist(1e-7, X_vx_ref, X_vx_c, L_ref, L_Earth_vx_BNS_c, 10, 10, events_count_vx_c, E_vx_BNS_c, 0.3, namelistc[2])

#%%

plt.scatter(E_MeV_BBH, X)
plt.show()

#%% Black Hole Merger 3x10^7 Msun
percent_count_BBH3_opt = plotlog(E_MeV_BBH, E_MeV_BBH, X_ve_ref, X, S_BBH, 814.8, L_ref, L_Earth_BBH_3_opt,'3x10^7 mass optimistic')
percent_count_BBH3_un = plotlog(E_MeV_BBH, E_MeV_BBH, X_ve_ref, X, S_BBH, 814.8, L_ref, L_Earth_BBH_3_un,'3x10^7 mass unoptimistic')

#Poisson Fluctuations
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X, L_Earth_BBH_3_opt, E_MeV_BBH, percent_count_BBH3_opt, '3x10^7 mass unoptimistic')
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X, L_Earth_BBH_3_un, E_MeV_BBH, percent_count_BBH3_un, '3x10^7 mass unoptimistic')

#Monte Carlo Sim
#energy_hist(1e-7, X_ve_ref, X_ve_i, L_ref, L_Earth_ve_BNS_i, 10, 10, events_count_ve_i, E_ve_BNS_i, 0.3, namelisti[0])
#energy_hist(1e-7, X_ve_ref, X_ve_c, L_ref, L_Earth_ve_BNS_c, 10, 10, events_count_ve_c, E_ve_BNS_c, 0.3, namelistc[0])

#%%Black Hole Merger 1x10^7 Msun
percent_count_BBH1_opt = plotlog(E_MeV_BBH, E_MeV_BBH, X_ve_ref, X, S_BBH, 814.8, L_ref, L_Earth_BBH_1_opt,'1x10^7 mass optimistic')
percent_count_BBH1_un = plotlog(E_MeV_BBH, E_MeV_BBH, X_ve_ref, X, S_BBH, 814.8, L_ref, L_Earth_BBH_1_un,'1x10^7 mass unoptimistic')

#Poisson Fluctuations 
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X, L_Earth_BBH_1_opt, E_MeV_BBH, percent_count_BBH1_opt, '1x10^7 mass unoptimistic')
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X, L_Earth_BBH_1_un, E_MeV_BBH, percent_count_BBH1_un, '1x10^7 mass unoptimistic')

#Monte Carlo Sim
#energy_hist(1e-7, X_ve_ref, X_ve_i, L_ref, L_Earth_ve_BNS_i, 10, 10, events_count_ve_i, E_ve_BNS_i, 0.3, namelisti[0])
#energy_hist(1e-7, X_ve_ref, X_ve_c, L_ref, L_Earth_ve_BNS_c, 10, 10, events_count_ve_c, E_ve_BNS_c, 0.3, namelistc[0])


#%% generate random luminosities, cross sections and distances

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

#%%
events_count_rand = plot(energy, E_1, X_ve_ref, X_E_1, S_d, 814.8, L_ref, Earth_Lum_1, 'Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 10e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X_vx_i, Earth_Lum_1, E_1, events_count_rand, 'Randomised Electron Neutrino')

#%% finding max distance event can be to still see neutrinos
Geo_trial = Torus(20)

S_trial = Geo_trial.Area()

#%%
#up to 94 for irro, 45.6 for coro
percent_count_ve_i = plot(energy, E_ve_BNS_i, X_ve_ref, Xve, S_trial, 814.8, L_ref, Geo_trial.Earth_Lum(L_ve_BNS_i * 1000),'Irrotational Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e51, 0.75e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])
percent_count_ve_c = plot(energy2, E_ve_BNS_c, X_ve_ref, Xve2, S_trial, 814.8, L_ref, Geo_trial.Earth_Lum(L_ve_BNS_c * 1000),'Corotational Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e51, 2e51], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

#%%
#up to 83 for irro, 40.5 for coro
events_count_vae_i = plot(energy, E_vae_BNS_i, X_vae_ref, Xvae, S_trial, 814.8, L_ref, Geo_trial.Earth_Lum(L_ve_BNS_i * 1000), 'Electron AntiNeutrino', t_BNS_i, [0.5, 30], [-0.1e52, 0.4e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])
events_count_vae_c = plot(energy2, E_vae_BNS_c, X_vae_ref, Xvae2, S_trial, 814.8, L_ref, Geo_trial.Earth_Lum(L_ve_BNS_c * 1000), 'Electron AntiNeutrino', t_BNS_c, [0.5, 30], [-0.1e52, 0.4e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

#%%
#up to 31.5 for irro, 28.5 for coro
events_count_vx_i = plot(energy, E_vx_BNS_i, X_vx_ref, Xvx, S_trial, 814.8, L_ref, Geo_trial.Earth_Lum(L_ve_BNS_i * 1000), 'Non-Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])
events_count_vx_c = plot(energy2, E_vx_BNS_c, X_vx_ref, Xvx2, S_trial, 814.8, L_ref, Geo_trial.Earth_Lum(L_ve_BNS_c * 1000), 'Non-Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 1e52], [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)])

#%% looking at the maximum distance the merger can be to still see neutrinos

distancelist = np.arange(1.0,100.0,5.0)

list_ve_i = []
list_ve_c = []
list_vae_i = []
list_vae_c = []
list_vx_i = []
list_vx_c = []

for i in range(distancelist.shape[0]):
    geo = Torus(distancelist[i])
    area = geo.Area()
    list_ve_i.append(distgraph(energy, E_ve_BNS_i, X_ve_ref, Xve, area, 814.8, L_ref, geo.Earth_Lum(L_ve_BNS_i * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_ve_c.append(distgraph(energy2, E_ve_BNS_c, X_ve_ref, Xve2, area, 814.8, L_ref, geo.Earth_Lum(L_ve_BNS_c * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))                                                                                                                                                                                            
    list_vae_i.append(distgraph(energy, E_vae_BNS_i, X_vae_ref, Xvae, area, 814.8, L_ref, geo.Earth_Lum(L_vae_BNS_i * 1000),[14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_vae_c.append(distgraph(energy2, E_vae_BNS_c, X_vae_ref, Xvae2, area, 814.8, L_ref, geo.Earth_Lum(L_vae_BNS_c * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_vx_i.append(distgraph(energy, E_vx_BNS_i, X_vx_ref, Xvx, area, 814.8, L_ref, geo.Earth_Lum(L_vx_BNS_i * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_vx_c.append(distgraph(energy2, E_vx_BNS_c, X_vx_ref, Xvx2, area, 814.8, L_ref, geo.Earth_Lum(L_vx_BNS_c * 1000),[14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))

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

#%%
angleslist = np.arange(0.0,np.pi/2,np.pi/40.0)

list_ve_i2 = []
list_ve_c2 = []
list_vae_i2 = []
list_vae_c2 = []
list_vx_i2 = []
list_vx_c2 = []

for i in range(angleslist.shape[0]):
    geo = Torus(50, angleslist[i])
    area = geo.Area()
    list_ve_i2.append(distgraph(energy, E_ve_BNS_i, X_ve_ref, Xve, area, 814.8, L_ref, geo.Earth_Lum(L_ve_BNS_i * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_ve_c2.append(distgraph(energy2, E_ve_BNS_c, X_ve_ref, Xve2, area, 814.8, L_ref, geo.Earth_Lum(L_ve_BNS_c * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))                                                                                                                                                                                            
    list_vae_i2.append(distgraph(energy, E_vae_BNS_i, X_vae_ref, Xvae, area, 814.8, L_ref, geo.Earth_Lum(L_vae_BNS_i * 1000),[14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_vae_c2.append(distgraph(energy2, E_vae_BNS_c, X_vae_ref, Xvae2, area, 814.8, L_ref, geo.Earth_Lum(L_vae_BNS_c * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_vx_i2.append(distgraph(energy, E_vx_BNS_i, X_vx_ref, Xvx, area, 814.8, L_ref, geo.Earth_Lum(L_vx_BNS_i * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
    list_vx_c2.append(distgraph(energy2, E_vx_BNS_c, X_vx_ref, Xvx2, area, 814.8, L_ref, geo.Earth_Lum(L_vx_BNS_c * 1000),[14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))

plt.plot(angleslist, list_ve_i2, marker = 'x', label = 'Irrotational ve', color = 'magenta')
plt.plot(angleslist, list_ve_c2, marker = 'x', label = 'Corotational ve', color = 'violet')
plt.plot(angleslist, list_vae_i2, marker = 'x', label = 'Irrotational vae', color = 'aqua')
plt.plot(angleslist, list_vae_c2, marker = 'x', label = 'Corotational vae', color = 'turquoise')
plt.plot(angleslist, list_vx_i2, marker = 'x', label = 'Irrotational vx', color = 'darkviolet')
plt.plot(angleslist, list_vx_c2, marker = 'x', label = 'Corotational vx', color = 'mediumpurple')
plt.legend()
plt.xlabel('Angles of Earth from Merger Plane (rad)')
plt.ylabel('Proportion of Events Above Threshold')
plt.show()

#%%

list_ve_i3 = []
list_ve_c3 = []
list_vae_i3 = []
list_vae_c3 = []
list_vx_i3 = []
list_vx_c3 = []

for i in range(distancelist.shape[0]):
    for j in range(angleslist.shape[0]):
         geo = Torus(distancelist[i], angleslist[j])
         area = geo.Area()
         list_ve_i3 = np.array( distgraph(energy, E_ve_BNS_i, X_ve_ref, Xve, area, 814.8, L_ref, geo.Earth_Lum(L_ve_BNS_i * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
         list_ve_c3.append(distgraph(energy2, E_ve_BNS_c, X_ve_ref, Xve2, area, 814.8, L_ref, geo.Earth_Lum(L_ve_BNS_c * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))                                                                                                                                                                                            
         list_vae_i3.append(distgraph(energy, E_vae_BNS_i, X_vae_ref, Xvae, area, 814.8, L_ref, geo.Earth_Lum(L_vae_BNS_i * 1000),[14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
         list_vae_c3.append(distgraph(energy2, E_vae_BNS_c, X_vae_ref, Xvae2, area, 814.8, L_ref, geo.Earth_Lum(L_vae_BNS_c * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
         list_vx_i3.append(distgraph(energy, E_vx_BNS_i, X_vx_ref, Xvx, area, 814.8, L_ref, geo.Earth_Lum(L_vx_BNS_i * 1000), [14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))
         list_vx_c3.append(distgraph(energy2, E_vx_BNS_c, X_vx_ref, Xvx2, area, 814.8, L_ref, geo.Earth_Lum(L_vx_BNS_c * 1000),[14, 20, 0.5e52, (-0.5e52 / 15), (0.5e52), (-0.5e52 / 12.5)]))


#%%
cf = conversion_factor(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10)
rate = np.array(X_ve_ref) * Earth_Lum_ref
g = cf * rate * 10
print(g)

plt.bar(9.5, g, color='turquoise', width = 0.1)
plt.show()

 #%%
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
#%%
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



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

t_SN, E_ve_SN, X_ve_SN, E_vae_SN, X_vae_SN, E_vx_SN, X_vx_SN, no_events, lum = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/energy_cross_sect.txt', 
                                                                skiprows = 1, unpack = True)
t_BNS_i, L_ve_BNS_i, E_ve_BNS_i, X_ve_BNS_i, L_vae_BNS_i, E_vae_BNS_i, X_vae_BNS_i, \
    L_vx_BNS_i, E_vx_BNS_i, X_vx_BNS_i = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BNS_data_irro.txt', skiprows = 1, unpack = True)

t_BNS_c, L_ve_BNS_c, E_ve_BNS_c, X_ve_BNS_c, L_vae_BNS_c, E_vae_BNS_c, X_vae_BNS_c, \
    L_vx_BNS_c, E_vx_BNS_c, X_vx_BNS_c = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/BNS_data_coro2.txt', skiprows = 1, unpack = True)

energy, Lve, Xve, Lvae, Xvae, Lvx, Xvx = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/energybins.txt', skiprows = 1, unpack = True)
energy2, Lve2, Xve2, Lvae2, Xvae2, Lvx2, Xvx2 = np.loadtxt('/Users/trinitystenhouse/Documents/University_MSci/2021-2/UROP/energybinscoro.txt', skiprows = 1, unpack = True)

#luminosity in ergs/s

Xve = Xve * 10 ** -38
Xvae = Xvae * 10 ** -38
Xvx = Xvx * 10 ** -38
 
Xve2 = Xve2 * 10 ** -38
Xvae2 = Xvae2 * 10 ** -38
Xvx2 = Xvx2 * 10 ** -38

Sphere_SN = Spherical(10) #input distance to event kpc
Sphere_BNS = Spherical(10)
S_SN = Sphere_SN.Area() #ref geometry
S_BNS = Sphere_BNS.Area()

L_ref_SN = max(Sphere_SN.Earth_Lum(lum))
L_Earth_ve_BNS_i = Sphere_BNS.Earth_Lum(L_ve_BNS_i * 1000)
L_Earth_vae_BNS_i = Sphere_BNS.Earth_Lum(L_vae_BNS_i * 1000)
L_Earth_vx_BNS_i = Sphere_BNS.Earth_Lum(L_vx_BNS_i * 1000)
L_Earth_ve_BNS_c = Sphere_BNS.Earth_Lum(L_ve_BNS_c * 1000)
L_Earth_vae_BNS_c = Sphere_BNS.Earth_Lum(L_vae_BNS_c * 1000)
L_Earth_vx_BNS_c = Sphere_BNS.Earth_Lum(L_vx_BNS_c * 1000)
#bin coro data

LEarthve = Sphere_BNS.Earth_Lum(Lve)
LEarthvae = Sphere_BNS.Earth_Lum(Lvae)
LEarthvx = Sphere_BNS.Earth_Lum(Lvx)

LEarthve2 = Sphere_BNS.Earth_Lum(Lve2)
LEarthvae2 = Sphere_BNS.Earth_Lum(Lvae2)
LEarthvx2 = Sphere_BNS.Earth_Lum(Lvx2)

Sphere_ref = Spherical(100) #using supernova at 100kpc
S_ref = Sphere_ref.Area()
L_ref = 0.5 * 10 ** 52
Earth_Lum_ref = max(Sphere_ref.Earth_Lum([L_ref]))
E_ref = 9.5
#10 interactions in 10s

X_ve_ref = 10 ** -5 * 10 ** -38
X_vae_ref = 13 ** -6 * 10 ** -38
X_vx_ref = 10 ** -6 * 10 ** -38

X_ve_i = X_ve_BNS_i * 10 ** -38
X_vae_i = X_vae_BNS_i * 10 ** -38
X_vx_i = X_vx_BNS_i * 10 ** -38 #all in cm^2
X_ve_c = X_ve_BNS_c * 10 ** -38
X_vae_c = X_vae_BNS_c * 10 ** -38
X_vx_c = X_vx_BNS_c * 10 ** -38 #all in cm^2

namelisti = ['Irrotational \n Electron Neutrino', 'Irrotational \n Electron AntiNeutrino', 'Irrotational \n Non Electron Neutrino']

namelistc = ['Corotational \n Electron Neutrino', 'Corotational \n Electron AntiNeutrino', 'Corotational \n Non Electron Neutrino']

#%% generate random luminosities, cross sections and distances

RP = random.randint(30,60)
RN = random.randint(1,100) * 10 ** RP
Space_Lum = ( L_ve_BNS_i / 10 ** 52 ) * RN #luminosity of an event at the event
RD = random.randint(1,20)
Sphere_Geo = Spherical(RD) #input distance to event
S_d = Sphere_Geo.Area() #ref geometry
Earth_Lum = Sphere_Geo.Earth_Lum(Space_Lum)
print(Earth_Lum)

X_E_1 = []
for i in range(0,len(X_ve_i)):
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
    
if np.median(Space_Lum) < np.median(lum):   
    print('Random Luminosities ', Space_Lum, 'ergs/s, less than')
else:
    print('Random Luminosities ', Space_Lum, 'ergs/s, more than')


#%%

Xvefull = np.concatenate((Xve, Xve2))
Xvaefull = np.concatenate((Xvae, Xvae2))
Xvxfull = np.concatenate((Xvx, Xvx2))
energyfull = np.concatenate((energy, energy2))

#%% events count and threshold plot for irrotational ideal
events_count_ve_i = plot(energy, E_ve_BNS_i, X_ve_ref, Xve, S_d, 100, L_ref, L_Earth_ve_BNS_i,'Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 4.5e52])

events_count_vae_i = plot(energy, E_vae_BNS_i, X_vae_ref, Xvae, S_d, 8, L_ref, L_Earth_vae_BNS_i, 'Electron AntiNeutrino', t_BNS_i, [0.5, 30], [-0.1e52, 8.5e52])

events_count_vx_i = plot(energy, E_vx_BNS_i, X_vx_ref, Xvx, S_d, 40, L_ref, L_Earth_vx_BNS_i, 'Non-Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 5.5e52])
#energy, E_v, X_ref, X_E, E, S_d, S_ref, L_ref, Earth_Lum, name, t

#%% Electron Neutrino
events_count_ve_i = plot(energy, E_ve_BNS_i, X_ve_ref, Xve, S_d, 1.9445e-17, L_ref, L_Earth_ve_BNS_i,'Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 1.5e52])
events_count_ve_c = plot(energy2, E_ve_BNS_c, X_ve_ref, Xve2, S_d, 1.9445e-17, L_ref, L_Earth_ve_BNS_c,'Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 4.5e52])

array_ve = plot2(energy, E_ve_BNS_i, E_ve_BNS_c, X_ve_ref, Xve, S_d, 1.9445e-17, L_ref, L_Earth_ve_BNS_i, L_Earth_ve_BNS_c, 'Electron Neutrino', t_BNS_i, t_BNS_c, [0.5, 30], [-0.1e52, 4.5e52])

#Poisson Fluctuations
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X_ve_i, L_Earth_ve_BNS_i, E_ve_BNS_i, namelisti[0])
events(np.array(X_ve_ref), np.array(Earth_Lum_ref), 10, 10, X_ve_c, L_Earth_ve_BNS_c, E_ve_BNS_c, namelistc[0])

#Monte Carlo Sim
#energy_hist(1e-7, X_ve_ref, X_ve_i, L_ref, L_Earth_ve_BNS_i, 10, 10, events_count_ve_i, E_ve_BNS_i, 0.3, namelisti[0])
#energy_hist(1e-7, X_ve_ref, X_ve_c, L_ref, L_Earth_ve_BNS_c, 10, 10, events_count_ve_c, E_ve_BNS_c, 0.3, namelistc[0])

#%% Electron AntiNeutrino
events_count_vae_i = plot(energy, E_vae_BNS_i, X_vae_ref, Xvae, S_d, 1.9445e-17, L_ref, L_Earth_vae_BNS_i, 'Electron AntiNeutrino', t_BNS_i, [0.5, 30], [-0.1e52, 8.5e52])
events_count_vae_c = plot(energy2, E_vae_BNS_c, X_vae_ref, Xvae2, S_d, 1.9445e-17, L_ref, L_Earth_vae_BNS_c, 'Electron AntiNeutrino', t_BNS_c, [0.5, 30], [-0.1e52, 8.5e52])

array_vae = plot2(energy, E_vae_BNS_i, E_vae_BNS_c, X_vae_ref, Xvae, S_d, 1.9445e-17, L_ref, L_Earth_vae_BNS_i, L_Earth_vae_BNS_c, 'Electron AntiNeutrino', t_BNS_i, t_BNS_c, [0.5, 30], [-0.1e52, 8.5e52])

#Poisson Fluctuations
events(np.array(X_vae_ref), np.array(Earth_Lum_ref), 10, 10, X_vae_i, L_Earth_vae_BNS_i, E_vae_BNS_i, namelisti[1])
events(np.array(X_vae_ref), np.array(Earth_Lum_ref), 10, 10, X_vae_c, L_Earth_vae_BNS_c, E_vae_BNS_c, namelistc[1])

#Monte Carlo Sim
#energy_hist(1e-7, X_vae_ref, X_vae_i, L_ref, L_Earth_vae_BNS_i, 10, 10, events_count_vae_i, E_vae_BNS_i, 0.27, namelisti[1])
#energy_hist(1e-7, X_vae_ref, X_vae_c, L_ref, L_Earth_vae_BNS_c, 10, 10, events_count_vae_c, E_vae_BNS_c, 0.3, namelistc[1])


#%% Non-Electron Neutrino 
events_count_vx_i = plot(energy, E_vx_BNS_i, X_vx_ref, Xvx, S_d, 1.9445e-17, L_ref, L_Earth_vx_BNS_i, 'Non-Electron Neutrino', t_BNS_i, [0.5, 30], [-0.1e52, 5.5e52])
events_count_vx_c = plot(energy2, E_vx_BNS_c, X_vx_ref, Xvx2, S_d, 1.9445e-17, L_ref, L_Earth_vx_BNS_c, 'Non-Electron Neutrino', t_BNS_c, [0.5, 30], [-0.1e52, 5.5e52])

array_vx = plot2(energy, E_vx_BNS_i, E_vx_BNS_c, X_vx_ref, Xvx, S_d, 1.9445e-17, L_ref, L_Earth_vx_BNS_i, L_Earth_vx_BNS_c, 'Non-Electron Neutrino', t_BNS_i, t_BNS_c, [0.5, 30], [-0.1e52, 5.5e52])

#Poisson Fluctuations
events(np.array(X_vx_ref), np.array(Earth_Lum_ref), 10, 10, X_vx_i, L_Earth_vx_BNS_i, E_vx_BNS_i, namelisti[2])
events(np.array(X_vx_ref), np.array(Earth_Lum_ref), 10, 10, X_vx_c, L_Earth_vx_BNS_c, E_vx_BNS_c, namelistc[2])

#Monte Carlo Sim
#energy_hist(1e-7, X_vx_ref, X_vx_i, L_ref, L_Earth_vx_BNS_i, 10, 10, events_count_vx_i, E_vx_BNS_i, 0.2, namelisti[2])
#energy_hist(1e-7, X_vx_ref, X_vx_c, L_ref, L_Earth_vx_BNS_c, 10, 10, events_count_vx_c, E_vx_BNS_c, 0.3, namelistc[2])

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

#%% make bins bigger
#%% github basic tutorial - version control 

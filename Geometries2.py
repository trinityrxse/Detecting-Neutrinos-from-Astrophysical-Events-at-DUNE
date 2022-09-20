#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 13:15:03 2022

@author: trinitystenhouse

Here is the module with the geometry classes, to specify the geometry of your event
"""

import numpy as np 

class Geometry:
    """
    allows us to call a geometry depending on event propagation
    """
    
class Spherical(Geometry):
    """ 
    represents geometry of an event modelled as propagation outward
    from a point source as a perfect sphere. 
    luminosities decrease according to inverse square law
    """
    def __str__(self):
         return f"Spherical Geometry of event, calculates luminosity \
             on Earth with respect to the geometry." 
         
    def __init__(self, r = 0):
        self._r = r #radius of sphere (distance to event) in kpc
    
    def Area(self):
        A = (4 * np.pi * (self._r ** 2))
        return A
        
    def Earth_Lum(self, Space_Lum):
        ELs = []
        for lum in Space_Lum:
            EL = lum / (self.Area())
            ELs.append(EL)
            
        return ELs

class Torus(Geometry):
    """ 
    represents geometry of an event modelled as propagation outward
    from a dual point source as two perfect spheres overlapping. 
    luminosities decrease with respect to angle of Earth wrt plane of event
    
    
    change this so it's ok again
    
    """
    def __str__(self):
        return f"Horn Torus Geometry of event, calculates luminosity \
            on Earth with respect to the geometry." 
    
    def __init__(self, D = 0, i = 0):
        self._D = D #distance to event in kpc
        self._i = i
    
    def Radius(self):
        r =  self._D / 2 * np.cos(self._i)
        return r
        
    def Area(self):
        A = ((4 * np.pi ** 2) * (self.Radius() ** 2))
        return A
    
    def Earth_Lum(self, Space_Lum):
        ELs = []
        for lum in Space_Lum:
            EL = lum / self.Area()
            ELs.append(EL)
            
        return ELs
 
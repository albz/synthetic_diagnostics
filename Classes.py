# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 00:04:05 2018

@author: Fil
"""
import numpy as np

class Vacuum:
    #define vacuum
    l = 0. #length [mm]
    def __init__(self, name):
        self.name = name
        
class Lens:    
    #define lens
    f=0. # focal length [mm]
    d=0. # aperture diameter [mm]
    L=0. # thickness [mm]
    def __init__(self, name):
        self.name = name
        
class Orifice:
    #define vacuum
    mask = []
    def __init__(self, name):
        self.name = name
        
class photon:
    #define photon for ray tracing
    x = []
    y = []
    z = []
    px= [] #momentum along x
    py= [] #momentum along y
    pz= [] #momentum along z
    wavelenght= 0.
    I=0.
#    def __init__(self, name):
#        self.name = name
    def __init__(self, x1, y1, z1, zdim, px1=0, py1=0, pz1=0, Intensity=1., wavelenght=0.) : #initialize the photon
        self.x=np.zeros(zdim)
        self.x[0]=x1
        self.y=np.zeros(zdim)
        self.y[0]=y1
        self.z=np.zeros(zdim)
        self.z[0]=z1
        self.px=np.zeros(zdim)
        self.px[0]=px1
        self.py=np.zeros(zdim)
        self.py[0]=py1
#        self.pz=np.zeros(zdim)
#        self.pz[0]=pz1
        self.I=Intensity
        self.wavelenght=wavelenght
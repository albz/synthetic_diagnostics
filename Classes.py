# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 00:04:05 2018

######################################################################
# Name:         Classes.py
# Author:		   F. Filippi
# Date:			2018-4-10
# Purpose:      classes for OpticPropagation.py
# Source:       python
#####################################################################
@author: Fil
"""
import numpy as np

class element:
    drift=0
    diopter=1
    thinlens=2
    mask=3
    grating=4
    
    def __init__(self, Type, L=0., R=0., f=0., n1=1., n2=1.,  gmm=1., maskType='none', D=0., matrix=np.ones(1)):
        self.type=Type
        if self.type == element.drift:
            self.length=L
            self.n=n1
        if self.type == element.diopter:
            self.Rcurvature=R
            self.n1=n1
            self.n2=n2
        if self.type == element.thinlens:
            self.focus=f
        if self.type == element.mask:
            self.masktype=maskType
            self.diameter=D
            self.matrix=matrix
        if self.type == element.grating:
            self.gmm=gmm
            
    def __str__(self):
        if self.type == element.drift:
            return 'drift'
        if self.type == element.diopter:
            return 'diopter'
        if self.type == element.thinlens:
            return 'thinlens'
        if self.type == element.mask:
            return 'mask'
        if self.type == element.grating:
            return 'grating'
 
        
class photon:
    #define photon for ray tracing
    x = []
    y = []
    z = []
    dx= [] #derivate along x
    dy= [] #derivate along y
    dz= [] #derivate along z
    wavelenght= 0.
    I=0.+0.j
#    def __init__(self, name):
#        self.name = name
    def __init__(self, x1, y1, z1, zdim, dx1=0, dy1=0, dz1=0, Intensity=1., wavelenght=0.) : #initialize the photon
        self.x=np.zeros(zdim)
        self.x[0]=x1
        self.y=np.zeros(zdim)
        self.y[0]=y1
        self.z=np.zeros(zdim)
        self.z[0]=z1
        self.dx=np.zeros(zdim)
        self.dx[0]=dx1
        self.dy=np.zeros(zdim)
        self.dy[0]=dy1
#        self.dz=np.zeros(zdim)
#        self.dz[0]=dz1
        self.I=Intensity
        self.wavelenght=wavelenght
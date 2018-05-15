#!/usr/bin/python
######################################################################
# Name:         OpticPropagation.py
# Author:		   F. Filippi
# Date:			2017-03-9
# Purpose:      optical geometry inteferometry from binary ALaDyn files
# Source:       python
#####################################################################
"""
Created on Thu Mar  9 15:10:30 2017

@author: Fil
"""
### loading shell commands
import os, os.path, sys
import numpy as np
import matplotlib.pyplot as plt

#sys.path.append(os.path.join(os.path.expanduser('~'),'Codes','Code_ALaDyn','tools-ALaDyn','pythons'))
#from read_ALaDyn_bin import *
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes','synthetic_diagnostics','tools','Functions'))
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes','synthetic_diagnostics','tools','Classes'))
from Functions import *
from Classes import *
#from sympy.physics.optics import RayTransferMatrix, ThinLens
#from sympy import Symbol, Matrix

"""
SETTINGS
"""
#NUMERICAL METHODS
analysis='Fourier' #'RayTracing'

#FLAGS
plotDensitySlicex=1;
plotDensitySlicey=0;
plotLaserSlice=0;
plotDephasingMap=0;
plotInterferogram=0;

#FUNDAMENTAL CONSTANTS
#c=299792458 #m/s
#e=1.60219E-19; #C
#me=9.1091E-31; #massa a riposo elettrone	kg 
#mp=1.6725E-27; #massa a riposo del protone	kg 
#kb=1.3806488E-23; #costante di Boltzmann J/K
#eps0=8.854187817e-12; #F/m

"""
INIZIALIZE PLASMA DENSITY
"""
#define the plasma density to cross
## - #
##LOAD FILE
#bin_path='C:/Users/Fil/Codes/synthetic_diagnostics'
#n_e,x,y,z = read_ALaDyn_bin(bin_path,'Edenout10.bin','grid')
#x=np.array(x)*1.e-6
#y=np.array(y)*1.e-6
#z=np.array(z)*1.e-6
#n_0=1E19 #CHIEDI ALBERTO NOME IN ALADYN
#n_e=n_e*n_0 #[cm^-3]
#
#ymax=y.size #high
#zmax=z.size #304 #probe laser propagation 
#xmax=x.size #1120 #main laser propagation
#Dx=0.2E-6; #m
#Dy=0.2E-6; #m
#Dz=0.2E-6; #m


#test density
ymax=300 #high
zmax=300 #304 #probe laser propagation 
xmax=300 #1120 #main laser propagation
Dx=1E-3; #m
Dy=1E-3; #m
Dz=1E-3; #m

x=np.linspace(0,Dx*xmax,xmax)
y=np.linspace(0,Dy*ymax,ymax)
z=np.linspace(0,Dz*zmax,zmax)

n=CreateTiltednMatrix(x,y,z,45,1.27)
#n_e=CreateTestDensity(xmax,ymax,zmax)    
n_0=1E19 #[cm^-3]
#n_e=n_e*n_0
    
print('Density acquired')

"""
INIZIALIZE PROBE LASER
"""
##define the probe profile
LambdaProbe=400e-9    #[m]
SpotxFWHM=1*Dx #200e-6      # m
SpotyFWHM=1*Dy #200e-6      # m
#LaserShape='gaussian' #or 'flattop'
#LaserSpot='gaussian'  #or 'flattop'
#Energy=1e-6           #in [J]
#PulseDuration=140e-15 #FWHM [s]
#Nx=400                #samples in x
#Ny=400                #samples in y
#Nz=1000               #samples in z
#xwindow=SpotxFWHM*2.123 #x window 5 sigma (5/2.35482=2.123)
#ywindow=SpotyFWHM*2.123 #y window 5 sigma
#zwindow=PulseDuration*c*2.123 #z window 5 sigma
#
#xprobe=np.linspace(-xwindow,xwindow, Nx) 
#yprobe=np.linspace(-ywindow,ywindow, Ny) 
#zprobe=np.linspace(-zwindow,zwindow, Nz)
#
##create the U(x,y,z) of the probe
##U [sqrt(W)/m]
#Upb=CreateLaserProfile(LaserShape,LaserSpot,LambdaProbe,Energy,PulseDuration,SpotxFWHM,SpotyFWHM, xprobe, yprobe, zprobe)

#Laser=CreateLaserphotonTest(x,y)

Laser=[photon]*int(SpotyFWHM/Dy)
for i in range(int(SpotxFWHM/Dx)):
    Laser[i]=[photon]*int(SpotxFWHM/Dx)
    for k in range(int(SpotyFWHM/Dy)):
        Laser[i][k]=photon(x[int(xmax/2+i)], y[int(ymax/2+k)], z[0], zmax)
        
        
print('Laser acquired')
    
"""
LOAD OPTICAL PATH
"""          
#how many lenses? how many BS? How to solve them (ray tracing or Fourier optics?) Where the plasma is?

"""
LASER PROPAGATION
"""
#yShift = np.fft.fftshift(y)
#fftyShift = np.fft.fft(yShift)
#ffty = np.fft.fftshift(fftyShift)

#convert plasma density in refrective index
#n_c=me*eps0*((2*np.pi*c)/e)**2/(LambdaProbe**2)*1e-6 #[cm^-3]
#n=np.sqrt(1-n_e/n_c)

#let the laser propagates inside 
zinit=z[0]
zend=z[-1]
Laser=RayTracingPropagator(x, y, z, Dx, Dy, Dz, zinit, zend, Laser, n)



     
"""
PLOTS
"""
#density
if plotDensitySlicex:
#    DensPlot = plt.figure(1) #, figsize=(3.25,3.0))
#    DPax  = plt.subplot(111)
    #fig.set_size_inches(2*3.25, 3.0, forward=True)
#    fig.set_size_inches(10.0, 5.0, forward=True)
    DPax=plt.plot(z, Laser[0][0].x, 'r')
    plt.imshow(n [:, int(ymax/2), :], extent=[z[0],z[-1],x[0],x[-1]], cmap=plt.cm.plasma)
    plt.xlabel('z')
    plt.ylabel('x')
    plt.colorbar()
    plt.show()
    
if plotDensitySlicey:
#    DensPlot = plt.figure(1) #, figsize=(3.25,3.0))
#    DPax  = plt.subplot(111)
    #fig.set_size_inches(2*3.25, 3.0, forward=True)
#    fig.set_size_inches(10.0, 5.0, forward=True)
    DPax=plt.plot(z, Laser[0][0].y, 'r')
    plt.imshow(n [int(xmax/2), :,  :], extent=[z[0],z[-1],y[0],y[-1]], cmap=plt.cm.plasma)
    plt.xlabel('z')
    plt.ylabel('y')
    plt.colorbar()
    plt.show()
    
#laser probe
if plotLaserSlice:
    LasPlot = plt.imshow(n [:, int(len(yprobe)/2), :], cmap=plt.cm.plasma)
    plt.colorbar()
    plt.xlabel('z')
    plt.ylabel('x')
    plt.show

#dephasing
if plotDephasingMap:
    DephPlot = plt.imshow(n [:, int(ymax/2), :], cmap=plt.cm.plasma)
    plt.colorbar
    plt.show
    
#interferogram
if plotInterferogram:
    IntegPlot = plt.imshow(interferogram.T, cmap=plt.cm.plasma)
    plt.show
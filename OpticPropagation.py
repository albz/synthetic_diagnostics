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
from scipy.constants import c

"""
SETTINGS
"""
#NUMERICAL METHODS
analysis='Fourier' #'RayTracing'

#FLAGS
plotDensitySlicex=0;
plotDensitySlicey=0;
plotLaserSlice=1;
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

xp=np.linspace(0,Dx*xmax,xmax)
yp=np.linspace(0,Dy*ymax,ymax)
zp=np.linspace(0,Dz*zmax,zmax)

#n=CreateTiltednMatrix(x,y,z,45,1.27)
n=CreateTiltednSlice(xp,yp,zp, 45, 100*Dx, 1.27)
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
LaserShape='CW'#'gaussian' #or 'flattop' or 'CW'
LaserSpot='flattop'#'gaussian'  #or 'flattop'
Energy=1e-6           #in [J]
PulseDuration=140e-15 #FWHM [s]
Nx=400                #samples in x
Ny=400                #samples in y
Nz=1000               #samples in z
xwindow=SpotxFWHM*2.123 #x window 5 sigma (5/2.35482=2.123)
ywindow=SpotyFWHM*2.123 #y window 5 sigma
zwindow=PulseDuration*c*2.123 #z window 5 sigma

x=np.linspace(-xwindow,xwindow, Nx) 
y=np.linspace(-ywindow,ywindow, Ny) 
z=np.linspace(-zwindow,zwindow, Nz)

#create the U(x,y,z) of the probe
#U [sqrt(W)/m]
Upb=CreateLaserProfile(LaserShape,LaserSpot,LambdaProbe,Energy,PulseDuration,SpotxFWHM,SpotyFWHM, x, y, z)

#Laser=CreateLaserphotonTest([SpotxFWHM/Dx], [SpotyFWHM/Dy], zp, zmax, lambdal=LambdaProbe)
#Laser=CreateLaserphotonTest([1,2,3], [1.5, 2.5, 3.5], zp, zmax, lambdal=LambdaProbe)
        
        
print('Laser acquired')
    
"""
LOAD OPTICAL PATH
"""
#how many lenses? how many BS? How to solve them (ray tracing or Fourier optics?) Where the plasma is?

"""
LASER PROPAGATION
"""

#convert plasma density in refrective index
#n_c=me*eps0*((2*np.pi*c)/e)**2/(LambdaProbe**2)*1e-6 #[cm^-3]
#n=np.sqrt(1-n_e/n_c)

#let the laser propagates inside 
zinit=zp[0]
zend=zp[-1]
#Laser=RayTracingPropagatorParax(xp, yp, zp, Dx, Dy, Dz, zinit, zend, Laser, n)
[x1, y1, Upbend]=Fresnel1StepPropagatorTF(x,y,Upb[:,:,0], LambdaProbe, 1)
#[x1, y1, upbend]=Fresnel2StepPropagatorTF(x,y,upb[:,:,0], LambdaProbe,0.1, 20)
#[x1, y1, upbend]=Fresnel2StepsPropagatorTF(x,y,upb[:,:,0], LambdaProbe,0.01, Dx) 



"""
PLOTS
"""
#density
if plotDensitySlicex:
#    DensPlot = plt.figure(1) #, figsize=(3.25,3.0))
#    DPax  = plt.subplot(111)
    #fig.set_size_inches(2*3.25, 3.0, forward=True)
#    fig.set_size_inches(10.0, 5.0, forward=True)
    DPax=plt.plot(zp, Laser[0][0].x, 'r')
    plt.imshow(n [:, int(ymax/2), :], extent=[zp[0],zp[-1],xp[0],xp[-1]], cmap=plt.cm.plasma)
    plt.xlabel('z')
    plt.ylabel('x')
    plt.colorbar()
    plt.show()
    
if plotDensitySlicey:
#    DensPlot = plt.figure(1) #, figsize=(3.25,3.0))
#    DPax  = plt.subplot(111)
    #fig.set_size_inches(2*3.25, 3.0, forward=True)
#    fig.set_size_inches(10.0, 5.0, forward=True)
    DPax=plt.plot(zp, Laser[0][0].y, 'r')
    plt.imshow(n [int(xmax/2), :,  :], extent=[zp[0],zp[-1],yp[0],yp[-1]], cmap=plt.cm.plasma)
    plt.xlabel('z')
    plt.ylabel('y')
    plt.colorbar()
    plt.show()
    
#laser probe
if plotLaserSlice:
    LasPlot = plt.imshow(Irradiance(Upbend), extent=[zp[0],zp[-1],yp[0],yp[-1]], cmap=plt.cm.plasma)
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
#!/usr/bin/python
######################################################################
# Name:         interferometry.py
# Author:		F. Filippi, A. Marocchino
# Date:			2017-07-13
# Purpose:      synthetic inteferometry from binary ALaDyn files
# Source:       python
#####################################################################

### loading shell commands
import os, os.path, sys
import math
import numpy as np
import matplotlib
import scipy.integrate as integrate
#matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import pylab as pyl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
### --- ###
sys.path.append(os.path.join(os.path.expanduser('~'),'Codes','Code_ALaDyn','tools-ALaDyn','pythons'))
from read_ALaDyn_bin import *
### --- ###

#FLAGS
plotDensitySlice=0;
plotDephasingMap=0;
plotInterferogram=1;


# - #
#LOAD FILE
bin_path='C:/Users/Fil/Codes/synthetic_diagnostics'
n,x,y,z = read_ALaDyn_bin(bin_path,'Edenout10.bin','grid')
n_0=1E19 #CHIEDI ALBERTO NOME IN ALADYN
n=n*n_0

#FUNDAMENTAL CONSTANTS
c=299792458 #[m/s]
e=1.60219E-19 # elementary charge [C]
me=9.1091E-31; # electron rest mass	[kg] 
eps0=8.854187817e-12; #[F/m]

#CCD 
#CCD_x=659
#CCD_z=494
recording_bit=16
MaxAmplitude=2**recording_bit

#PROBE BEAM PROPERTIES
LambdaProbe=400e-9 #[m]
PercentSignal=100 #[%] over MaxAmplitude
PercentNoise=0 #[%] over MaxAmplitude
SignalAmpl=MaxAmplitude*PercentSignal*1E-2

#INTERFEROMETER SETUP
LSpatial_z=math.inf  #pixels per fringe along z axis
LSpatial_x=40 #pixels per fringe along x axis





###CALCULATIONS###
n_c=me*eps0*((2*math.pi*c)/e)**2/(LambdaProbe**2)*1e-6 #[cm^-3]

Dphi=(2*math.pi/LambdaProbe)*integrate.trapz(1-np.sqrt(1-n/n_c), np.array(y,dtype=float)*1.e-6, axis=1) #dephasing map (integrated along y axis)

interferogram=np.zeros(Dphi.shape)


for i in range(1,len(x)):
    for j in range(1,len(z)):
        interferogram[i,j]=SignalAmpl*math.cos(2*math.pi*i/LSpatial_x + 2*math.pi*j/LSpatial_z + Dphi[i,j])


#PLOT
#density
if plotDensitySlice:
    DensPlot = plt.imshow(n [:, len(y)/2, :])
    plt.colorbar
    plt.show

#dephasing
if plotDephasingMap:
    DephPlot = plt.imshow(n [:, len(y)/2, :])
    plt.colorbar
    plt.show
    
#interferogram
if plotInterferogram:
    IntegPlot = plt.imshow(interferogram.T)
    plt.show
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 13:40:06 2018

@author: Fil
"""

def CreateLaserProfile(LaserShape,LaserSpot,LambdaProbe,Energy,PulseDuration,SpotxFWHM,SpotyFWHM, x,y,z):
    #x,y,z in [m]
    # 
    
    c=299792458 #m/s
    #temporal profile
    TemporalProfile=np.zeros(len(z))
    if LaserShape=='gaussian':
        sigmaz=c*PulseDuration*2*np.sqrt(2*np.log(2)) #[m]
        az=1/(sigmaz * np.sqrt(2*np.pi))
        muz=0;
        
        TemporalProfile=az*np.exp(-1/2*((z-muz)/sigmaz)**2)
    if LaserShape=='flattop':
        TemporalProfile=[1/(c*PulseDuration) if c else 0 for(c) in [((z>-(c*PulseDuration)) & (z<(c*PulseDuration))) for z in z]]
        
        
    #transverse profile
    SpotProfile=np.zeros(len(x),len(y))
    if LaserSpot=='gaussian':
        sigmax=SpotxFWHM*2*np.sqrt(2*np.log(2)) #[m]
        ax=1/(sigmax * np.sqrt(2*np.pi))
        mux=0;
        
        sigmay=SpotyFWHM*2*np.sqrt(2*np.log(2)) #[m]
        ay=1/(sigmay * np.sqrt(2*np.pi))
        muy=0;
        
        SpotProfile=ax*ay*np.exp(-1/2*(((x-mux)/sigmax)**2+((y-muy)/sigmay)**2))
    if LaserSpot=='flattop':
        #TO BE IMPLEMENTED
        SpotProfile=1/SpotxFWHM;
        
        
    LaserProfile
    
    return LaserProfile
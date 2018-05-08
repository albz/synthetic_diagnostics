# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:58:12 2018

@author: Fil
"""

def CreateTiltednMatrix(x,y,z, angletox, nin):
    #x,y,z coordinates of the matrix n
    #z is the light propagation axis
    #angletox angle in degree with the x axis (perpendicular to z)
    #n refractive index of the medium
    n=np.ones([len(x),len(y),len(z)])
    
    centerx=np.abs(x[0]-x[-1])/2
    centerz=np.abs(z[0]-z[-1])/2
#    print(centerx)
#    print(centerz)
            
    for i in range(len(x)):
        for j in range(len(z)):
            if (z[j]-centerz)*np.sin(np.radians(90-angletox))>(x[i]-centerx)*np.cos(np.radians(90-angletox)):
#            if x[i]<(z[j]-z[1])*np.tan(np.radians(angletox)):
                n[i,:,j]=nin
    
    return n
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 13:38:37 2018

@author: Fil
"""

def CreateTestDensity(x,y,z):
    n=np.zeros([x,y,z])
    #xquarter=int(x/4)
    #yquarter=int(y/4)
    #n[xquarter:x-xquarter, yquarter:y-yquarter, :]=np.ones([xquarter*2,yquarter*2,z])
    for i in range(1,x):
        for j in range(1,y):
            for k in range(1,z):
                if k<i*(z/x):
                    n[i,j,k]=1.
                     
    return n

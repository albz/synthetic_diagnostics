# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 18:30:37 2017

@author: Fil
"""
import numpy as np
from scipy import interpolate
import scipy.integrate as integrate
from scipy.integrate import odeint
import matplotlib.pyplot as plt



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

def CreateTiltednMatrix(x,y,z, angletox, nin):
    #x,y,z coordinates of the matrix n
    #z is the light propagation axis
    #angletox angle in degree with the x axis (perpendicular to z)
    #n refractive index of the medium
    n=np.ones([len(x),len(y),len(z)])
    
    centerx=np.abs(x[0]-x[-1])/2
    centerz=np.abs(z[0]-z[-1])/2
            
    for i in range(len(x)):
        for j in range(len(z)):
            if (z[j]-centerz)*np.sin(np.radians(90-angletox))>(x[i]-centerx)*np.cos(np.radians(90-angletox)):
                n[i,:,j]=nin    
    return n


def CreateLaserProfile(LaserShape,LaserSpot,LambdaProbe,Energy,PulseDuration,SpotxFWHM,SpotyFWHM, x,y,z):
    #x,y,z in [m]
    # LaserProfile [sqrt(W)/m] (it is the Upb=[Upbx, Upby, Upbz])
    
#    c=299792458 #m/s
    #temporal profile
    TemporalProfile=np.zeros(len(z))
    if LaserShape=='gaussian':
        sigmaz=c*PulseDuration*2*np.sqrt(2*np.log(2)) #[m]
        az=1/(sigmaz * np.sqrt(2*np.pi))
        muz=0;
        
        TemporalProfile=az*np.exp(-1/2*((z-muz)/sigmaz)**2)
    if LaserShape=='flattop':
        TemporalProfile=[1/(c*PulseDuration) if c else 0 for(c) in [((z>-(c*PulseDuration)) & (z<(c*PulseDuration))) for z in z]]
    
    Upbz=TemporalProfile    
        
    #transverse profile
    SpotProfile=np.zeros((len(x),len(y)))
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
        
        
#    LaserProfile
    
    #LambdaProbe,
    #normalize to energy
    TempInt=np.trapz(np.trapz(np.trapz(LaserProfile**2, x, axis=0), y, axis=1), z/c, axis=2)
    LaserProfile=LaserProfile*Energy/TempInt
    
    return Upbz


"""
FFT
"""
def fft(f, delta):
    #f function to transform 1D
    #delta step of the sampling (assumed equispaced)
    F=np.fft.fftshift(np.fft.fft((f)))*delta
    return F

def ifft(F, deltaF):
    #F function to inverse transform 1D
    #deltaF step of the sampling (assumed equispaced)
    f=np.fft.ifftshift(np.fft.ifft((F)))*deltaF
    return f

def fft2(f, delta):
    #f function to transform 2D
    #delta step of the sampling (assumed equispaced)
    F=np.fft.fftshift(np.fft.fft2((f)))*delta**2
    return F

def ifft2(F, deltaF):
    #F function to inverse transform 2D
    #deltaF step of the sampling (assumed equispaced)
    f=np.fft.ifftshift(np.fft.ifft2((F)))*deltaF**2
    return f


"""
DERIVATIVES
"""
def gradient_ft(g, delta):
    #gradient operation based on the fft
    #g is 2D function
    N=np.size(g,0) #samples per size
    F=1/(N*delta)
    fX=np.range(-N/2,N/2-1)*F
    fX, fY=np.meshgrid(fX)
    gx=ifft2(-2j*np.pi*fX*fft2(g,delta),F)
    gx=ifft2(-2j*np.pi*fY*fft2(g,delta),F)
    return [gx, gy]
       

"""
SOLVERS
"""
def FresnelPropagatorTF(x,y,Uin, Htf, wavelength, z) :
    #Uin source plane field
    #Htf transformed transfer function
    #wavelength is given in [m] 
    #z distance of propagation in [m]
    from scipy.constants import pi
    Uin=fft2(uin, delta):
    k=2*pi/wavelength #wavenumber
    
    
    return Uout

def TFVacuum(z) :
    #z propagation distance
    Uin
    
    
    return Uvac

def FraunhoferPropagatorTF(Uin, Htf, wavelength) :
    #Uin source plane field
    #Htf transformed transfer function
    Uin
    
    
    return Uout

def ConvertPlasma2n(n_e, LambdaProbe):
    #convert plasma density in refrective index
    #LambdaProbe photon wavelength in [m]
    from scipy.constants import c,e, m_e, epsilon_0, pi
#    c=299792458 #m/s
#    e=1.60219E-19; #C
#    m_e=9.1091E-31; #massa a riposo elettrone	kg 
#    mp=1.6725E-27; #massa a riposo del protone	kg 
#    kb=1.3806488E-23; #costante di Boltzmann J/K
#    eps0=8.854187817e-12; #F/m
    
    n_c=m_e*epsilon_0*((2*pi*c)/e)**2/(LambdaProbe**2)*1e-6 #[cm^-3]
    n=np.sqrt(1-n_e/n_c)
    return n


def F(yin, z, dndx, dndy, dndz, n): 
    # Return derivatives for second-order ODE y'' = - a(z) y' + b(z)
    x=[0, 0]
    y=[0, 0]
    x[0]=yin[0] #x(t)
    x[1]=yin[1] #first derivative of x(t).
    y[0]=yin[2] #y(t)
    y[1]=yin[3] #first derivative of y(t).
    
    dyin = [0, 0, 0, 0]    # Create a list to store derivatives.
    dyin[0] = yin[1]    # Store first derivative of x(t).
    dyin[1] = -dndz([x[0],y[0],z])/n([x[0],y[0],z])*x[1]+dndx([x[0],y[0],z])/n([x[0],y[0],z])    # Store second derivative of x(t).
    dyin[2] = y[1]    # Store first derivative of y(t).
    dyin[3] = -dndz([x[0],y[0],z])/n([x[0],y[0],z])*y[1]+dndy([x[0],y[0],z])/n([x[0],y[0],z])    # Store second derivative of y(t).
    return dyin


def Fx(xin, z, a, b): 
    # Return derivatives for second-order ODE y'' = - a(z) y' + b(z)
    x=[0, 0]
    x[0]=xin[0] #x(t)
    x[1]=xin[1] #first derivative of x(t)
    
    dxin = [0, 0]    # Create a list to store derivatives.
    dxin[0] = xin[1]    # Store first derivative of x(t).
    dxin[1] = -1*a*x[1]+b    # Store second derivative of x(t).
    return dxin


def Fy(yin, z, a, b): 
    # Return derivatives for second-order ODE y'' = - a(z) y' + b(z)
    y=[0, 0]
    y[0]=yin[0] #y(t)
    y[1]=yin[1] #first derivative of y(t)
    
    dyin = [0, 0]    # Create a list to store derivatives.
    dyin[0] = yin[1]    # Store first derivative of y(t).
    dyin[1] = -1*a*y[1]+b    # Store second derivative of y(t).
    return dyin



def RayTracingPropagator(x, y, z, Dx, Dy, Dz, zinit, zend, RayMatrix, n) :
    #x, y, z the coordinates of the the matrix n along which the light is propagating
    #NB z is the light propagation axis
    #zinit is the integration start
    #zend is the integration end
    #RayMatrix is the matrix of ray elements of class 'photon' described in Classes.py
    #n is the refractive index matrix  
    #THIS VERSION DERIVATES BEFORE PROPAGATING
    
    from scipy.interpolate import RegularGridInterpolator as rgi
    
    #evaluates the derivatives of n
    dndx=np.diff(n,axis=0)/Dx
    dndxi = rgi((x[0:-1]+Dx/2, y, z), dndx, bounds_error=False, fill_value=None)
    
    dndy=np.diff(n, axis=1)/Dy
    dndyi = rgi((x, y[0:-1]+Dy/2, z), dndy, bounds_error=False, fill_value=None)
    
    dndz=np.diff(n, axis=2)/Dz
    dndzi = rgi((x, y, z[0:-1]+Dz/2), dndz, bounds_error=False, fill_value=None)   
    
    ni = rgi((x, y, z), n, method='linear', bounds_error=False, fill_value=1)
    
    init=np.argmin(np.abs(z-zinit))
    stop=np.argmin(np.abs(z-zend))
    
    for i in range(len(RayMatrix)):
        for j in range(len(RayMatrix[0])):
            print('Ray Tracing')
            print(((len(RayMatrix)-i)/len(RayMatrix))*100, '%')
            print('\n')
            
            #it analyse separatly x and y
            for zcount in range(init,stop):
                n_local=ni([RayMatrix[i][j].x[zcount],RayMatrix[i][j].y[zcount],z[zcount]])
                x0=[RayMatrix[i][j].x[zcount],  RayMatrix[i][j].px[zcount]/n_local] #initial conditions x and x'
                y0=[RayMatrix[i][j].y[zcount],  RayMatrix[i][j].py[zcount]/n_local] #initial conditions y and y'
                
                dndx_local=dndxi([x0[0],y0[0],z[zcount]])
                dndy_local=dndyi([x0[0],y0[0],z[zcount]])
                dndz_local=dndzi([x0[0],y0[0],z[zcount]])                
                
                
                xsol = odeint(Fx, x0, z[range(zcount, zcount+2)], args=(dndz_local/n_local, -dndx_local/n_local), mxstep=5000000)#,mxstep=5000000
                ysol = odeint(Fy, y0, z[range(zcount, zcount+2)], args=(dndz_local/n_local, -dndy_local/n_local), mxstep=5000000)
                
#                if dndx_local!=0: 
#                    print(dndx_local, dndy_local, dndz_local, n_local, xsol[1,:], ysol[1,:], z[zcount])
                
                RayMatrix[i][j].x[zcount+1]=xsol[1,0]
                RayMatrix[i][j].y[zcount+1]=ysol[1,0]  
                RayMatrix[i][j].z[zcount+1]=z[zcount+1]
                RayMatrix[i][j].px[zcount+1]=xsol[1,1]*ni([xsol[1,0],ysol[1,0], z[zcount+1]])
                RayMatrix[i][j].py[zcount+1]=ysol[1,1]*ni([xsol[1,0],ysol[1,0], z[zcount+1]])            
                
    return RayMatrix

def RayTracingPropagator2V(x, y, z, Dx, Dy, Dz, zinit, zend, RayMatrix, n) :
    #x, y, z the coordinates of the the matrix n along which the light is propagating
    #NB z is the light propagation axis
    #zinit is the integration start
    #zend is the integration end
    #RayMatrix is the matrix of ray elements of class 'photon' described in Classes.py
    #n is the refractive index matrix  
    #THIS VERSION USE THE SPLINE INTERPOLATION TO EVALUATE THE DERIVATIVES
    
    from scipy.interpolate import RectBivariateSpline as rbs
        
    init=np.argmin(np.abs(z-zinit))
    stop=np.argmin(np.abs(z-zend))
    
    for i in range(len(RayMatrix)):
        for j in range(len(RayMatrix[0])):
            print('Ray Tracing')
            print(((len(RayMatrix)-i)/len(RayMatrix))*100, '%')
            print('\n')

            for zcount in range(init,stop): 
                ycount=np.amin(RayMatrix[i][j].y[zcount]-y)
                ni=rbs(x,y, n [:, :, zcount])
                niz=rbs(x,z, n [:, ycount, :])
                x0=[RayMatrix[i][j].x[zcount],  RayMatrix[i][j].px[zcount]/ni(RayMatrix[i][j].x[zcount],RayMatrix[i][j].y[zcount])] #initial conditions x and x'
                y0=[RayMatrix[i][j].y[zcount],  RayMatrix[i][j].py[zcount]/ni(RayMatrix[i][j].x[zcount],RayMatrix[i][j].y[zcount])] #initial conditions y and y'
                
                dndx_local=ni(x0[0],y0[0],dx=1)
                dndy_local=ni(x0[0],y0[0],dy=1)
                dndz_local=niz(x0[0],z[zcount],dy=1) 
                n_local=ni(x0[0],y0[0])
                
                
                xsol = odeint(Fx, x0, z[range(zcount, zcount+2)], args=(dndz_local/n_local, -dndx_local/n_local), mxstep=5000000)#,mxstep=5000000
                ysol = odeint(Fy, y0, z[range(zcount, zcount+2)], args=(dndz_local/n_local, -dndy_local/n_local), mxstep=5000000)
                
#                if dndx_local!=0: 
#                    print(dndx_local, dndy_local, dndz_local, n_local, xsol[1,:], ysol[1,:], z[zcount])
                
                RayMatrix[i][j].x[zcount+1]=xsol[1,0]
                RayMatrix[i][j].y[zcount+1]=ysol[1,0]  
                RayMatrix[i][j].z[zcount+1]=z[zcount+1]
                RayMatrix[i][j].px[zcount+1]=xsol[1,1]*n_local
                RayMatrix[i][j].py[zcount+1]=ysol[1,1]*n_local    
                         
                
    return RayMatrix
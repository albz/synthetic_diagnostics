# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 18:30:37 2017

######################################################################
# Name:         Functions.py
# Author:		   F. Filippi
# Date:			2017-3-9
# Purpose:      classes for OpticPropagation.py
# Source:       python
#####################################################################

@author: Fil
"""
from Classes import *



def CreateTestDensity(x,y,z):
    import numpy as np
    
    n=np.zeros([x,y,z])
    #xquarter=int(x/4)
    #yquarter=int(y/4)
    for i in range(1,x):
        for j in range(1,y):
            for k in range(1,z):
                if k<i*(z/x):
                    n[i,j,k]=1.
                   
    return n


def CreateTiltednMatrix(x,y,z, angletox, nin):
    """
    This function create a matrix tilted with an angle given in degrees
    x,y,z coordinates of the matrix n
    z is the light propagation axis
    angletox angle in degree with the x axis (perpendicular to z)
    n refractive index of the medium
    """
    n=np.ones([len(x),len(y),len(z)])
    
    centerx=np.abs(x[0]-x[-1])/2
    centerz=np.abs(z[0]-z[-1])/2
            
    for i in range(len(x)):
        for j in range(len(z)):
            if (z[j]-centerz)*np.sin(np.radians(90-angletox))>(x[i]-centerx)*np.cos(np.radians(90-angletox)):
                n[i,:,j]=nin    
    return n


def CreateTiltednSlice(x,y,z, angletox, thickness, nin):
    """
    x,y,z coordinates of the matrix n
    z is the light propagation axis
    angletox angle in degree with the x axis (perpendicular to z)
    thickness of the slice
    n refractive index of the medium
    """
    import numpy as np
    
    n=np.ones([len(x),len(y),len(z)])
    
    centerx=np.abs(x[0]-x[-1])/2
    centerz=np.abs(z[0]-z[-1])/2
    centerx2=thickness*np.cos(np.radians(90-angletox))+centerx
    centerz2=thickness*np.sin(np.radians(90-angletox))+centerz
#    print(centerx2,centerx,centerz2,centerz)
            
    for i in range(len(x)):
        for j in range(len(z)):
#            if (z[j]-centerz2)*np.sin(np.radians(90-angletox))>(x[i]-centerx2)*np.cos(np.radians(90-angletox)):
            if (z[j])*np.sin(np.radians(90-angletox))>(x[i])*np.cos(np.radians(90-angletox)):
                if (z[j]-centerz2)*np.sin(np.radians(90-angletox))<(x[i]-(centerx2-thickness))*np.cos(np.radians(90-angletox)):
                    n[i,:,j]=nin; 
    return n


def CreateLaserProfile(LaserShape,LaserSpot,LambdaProbe,Energy,PulseDuration,SpotxFWHM,SpotyFWHM, x,y,z):
    """
    This function creates a 3D matrix 
    
    LaserProfile [sqrt(W)/m] (it is the fasor matrix Upb=[Upbxy, Upbz])
    
    LaserShape    can be 'gaussian', 'flattop' or continuous 'CW'
    LaserSpot    can be 'gaussian', 'flattop'
    LambdaProbe [m]
    Energy      [J]
    PulseDuration    FWHM [s]
    SpotxFWHM   [m]
    SpotyFWHM   [m]
    x,y,z in [m]
    """
    import numpy as np
    from scipy.constants import c
    
    if (LambdaProbe>abs(z[1]-z[0])) & (LaserShape!='CW'):
        print('WARNING! Temporal resolution is finer than the wavelength')
    
    #temporal profile
    TemporalProfile=np.zeros(len(z))
    if LaserShape=='gaussian':
        sigmaz=c*PulseDuration*2*np.sqrt(2*np.log(2)) #[m]
        az=1/(sigmaz * np.sqrt(2*np.pi))
        muz=0;
        
        TemporalProfile=az*np.exp(-1/2*((z-muz)/sigmaz)**2)
    if LaserShape=='flattop':
        TemporalProfile=[1/(c*PulseDuration) if c else 0 for(c) in [((z>-(c*PulseDuration)) & (z<(c*PulseDuration))) for z in z]]
    
    if LaserShape=='CW':
        TemporalProfile=[1.]
        
        
    #transverse profile
    SpotProfile=np.zeros((len(x),len(y)))
    if LaserSpot=='gaussian':
        sigmax=SpotxFWHM*2*np.sqrt(2*np.log(2)) #[m]
        ax=1/(sigmax * np.sqrt(2*np.pi))
        mux=0;        
        sigmay=SpotyFWHM*2*np.sqrt(2*np.log(2)) #[m]
        ay=1/(sigmay * np.sqrt(2*np.pi))
        muy=0;        
        for i in range(np.size(x)):
            SpotProfile[i,:]=ax*ay*np.exp(-1/2*(((x[i]-mux)/sigmax)**2+((y-muy)/sigmay)**2))
            
    if LaserSpot=='flattop':
        Xprofile=[1/(SpotxFWHM) if c else 0 for(c) in [((x>-(SpotxFWHM)) & (x<(SpotxFWHM))) for x in x]]
        Yprofile=[1/(SpotyFWHM) if c else 0 for(c) in [((y>-(SpotyFWHM)) & (y<(SpotyFWHM))) for y in y]]          
        for i in range(np.size(x)):
            SpotProfile[i,:]=[Xprofile[i]*x for x in Yprofile]
        
#    create LaserProfile
    LaserProfile=np.zeros((len(x),len(y),np.size(TemporalProfile)))
    for i in range(np.size(TemporalProfile)):
            LaserProfile[:,:,i]=[TemporalProfile[i]*x for x in SpotProfile]
    
    #Set the energy of the laser pulse
    LaserProfile=LaserProfile*Energy/LaserEnergy(LaserProfile, x,y,z)
    
    return LaserProfile

def CreateLaserphotonTest(x,y, z, zend, lambdal=0.):
    """
    It creates a 3D matrix of photon element
    """
    import numpy as np
    
    Laser=[photon]*np.size(x)
    for i in range(int(np.size(x))-1):
        Laser[i]=[photon]*int(np.size(y))
        for k in range(int(np.size(y))-1):
            print(i,k)
            Laser[i][k]=photon(x[int(np.size(x)/2+i)], y[int(np.size(y)/2+k)], z[0], zend, wavelenght=lambdal)
    return Laser

"""
FFT
"""
def fft(f, delta):
    """
    f function to transform 1D
    delta step of the sampling (assumed equispaced)
    """
    import numpy as np
    
    F=np.fft.fftshift(np.fft.fft((f)))*delta
    return F

def ifft(F, deltaF):
    """
    F function to inverse transform 1D
    deltaF step of the sampling (assumed equispaced)
    """
    import numpy as np
    
    f=np.fft.ifftshift(np.fft.ifft((F)))*deltaF
    return f

def fft2(f, delta):
    #f function to transform 2D
    #delta step of the sampling (assumed equispaced)
    import numpy as np
    
    F=np.fft.fftshift(np.fft.fft2((f)))*delta**2
    return F

def ifft2(F, deltaF):
    #F function to inverse transform 2D
    #deltaF step of the sampling (assumed equispaced)
    import numpy as np
    
    f=np.fft.ifftshift(np.fft.ifft2((F)))*deltaF**2
    return f


"""
DERIVATIVES
"""
def gradient_ft(g, delta):
    """
    gradient operation based on the fft
    g is 2D matrix of float value
    delta is the sampling step (assumed equispaced)
    """
    import numpy as np
    
    N=np.size(g,0) #samples per size
    F=1/(N*delta)
    fX=np.linspace(-N/2,N/2-1,N)*F
    fX, fY=np.meshgrid(fX, fX)
    gx=ifft2(-2j*np.pi*fX*fft2(g,delta),F)
    gy=ifft2(-2j*np.pi*fY*fft2(g,delta),F)
    return [gx, gy]
       

"""
SOLVERS
"""
#FOURIER OPTICS
def Fresnel1StepPropagatorTF(x,y,Uin, wavelength, Dz, n=1., printdelta='False') :
    """
    Fresnel one step propagator for monochromatic wave in homogenous medium
    uin source plane field [sqrt(W)/m]
    wavelength is given in [m] 
    Dz distance of propagation in [m]
    n is the refractive index sqrt(epsilon/epsilon_0), default 1., [adim]
    printdelta if True it prints the spacing grid of the final spot
    """
    import numpy as np
    from scipy.constants import pi    
    
    N = np.size(x)
    d1x=np.abs(x[-1]-x[0])/N #x grid spacing of the initial spot
    d2x=wavelength*Dz/(N*d1x) #x grid spacing of the final spot
    M = np.size(y)
    d1y=np.abs(y[-1]-y[0])/M
    d2y=wavelength*Dz/(M*d1y)
    if printdelta=='True':
        print('dxobj=',d2x,'dyobj=',d2y)
    
    k = n*2*pi/wavelength;    # optical wavevector    
    x1, y1 = np.meshgrid(x,y) # source-plane coordinates
    x2, y2 = np.meshgrid(np.linspace(-N/2,N/2-1,N)*d2x, np.linspace(-M/2,M/2-1,M)*d2y)  # observation-plane coordinates
    # evaluate the Fresnel-Kirchhoff integral
    Uout = 1/(1j*wavelength*Dz) * np.exp(1j*k/(2*Dz)*(x2**2 + y2**2)) * fft2(Uin * np.exp(1j*k/(2*Dz)*(x1**2 + y1**2)), d1x)
    
    return [x2, y2, Uout]

def Fresnel2StepPropagatorTF(x,y,Uin, wavelength, Dz, d2, n=1., printdelta='False') :
    """
    Fresnel two steps propagator for monochromatic wave in homogenous medium
    uin source plane field [sqrt(W)/m]
    wavelength is given in [m] 
    d2x, d2y final grid spacing 
    Dz distance of propagation in [m]
    n is the refractive index sqrt(epsilon/epsilon_0), default 1., [adim]
    printdelta if True it prints the spacing grid of the initial and final spot
    """
    import numpy as np
    from scipy.constants import pi    
    
    N = np.size(x) # number of grid points
    d1x=np.abs(x[-1]-x[0])/N #x grid spacing of the initial spot
    M = np.size(y)
    d1y=np.abs(y[-1]-y[0])/M          
              
    # source-plane coordinates
    k = n*2*pi/wavelength;    # optical wavevector    
    x1, y1 = np.meshgrid(x,y) # source-plane coordinates
    
    # magnification
    m = d2/d1x;
    # intermediate plane
    Dz1 = Dz / (1 - m); # propagation distance 
    d1a = wavelength * abs(Dz1) / (N * d1x);
    # coordinates
    x1a, y1a = np.meshgrid(x*d1a/d1x,y*d1a/d1y) # intermidiate plane coordinates
    #evaluate the Fresnel-Kirchhoff integral
    Uitm = 1/(1j*wavelength*Dz1) * np.exp(1j*k/(2*Dz1) * (x1a**2+y1a**2))* fft2(Uin*np.exp(1j*k/(2*Dz1)*(x1**2 + y1**2)), d1x)
    # observation plane
    Dz2 = Dz - Dz1; # propagation distance # coordinates
    x2, y2 = np.meshgrid(x*d2/d1x,y*d2/d1y) # source-plane coordinates
    # evaluate the Fresnel diffraction integral
    Uout = 1/(1j*wavelength*Dz2) * np.exp(1j*k/(2*Dz2) * (x2**2+y2**2)) * fft2(Uitm*np.exp(1j*k/(2*Dz2) * (x1a**2 + y1a**2)), d1a)
    
    if printdelta:
        print('dxobj=',d2,'dyobj=',d2)
        
    return [x2, y2, Uout]


def TFVacuum(z) :
    #z propagation distance
    Uin
    
    
    return Uvac

def FraunhoferPropagatorTF(Uin, Htf, wavelength) :
    #Uin source plane field
    #Htf transformed transfer function
    Uin
    
    
    return Uout


#GEOMETRICAL OPTIC
def OpticalMatrices(OptElement) :
    """
    This function return the optical matrix 2X2 relative to the input of class
    element          
    """
    import numpy as np
    
    if OptElement.type == element.drift:
            return [[1., (OptElement.L)],[0., 1.]] #DA FARE
    if OptElement.type == element.diopter:
            return [[1., 0.],[0., 1.]] #DA FARE
    if OptElement.type == element.thinlens:
            return [[1., 0.],[-1./OptElement.focus, 1.]] #DA FARE
    if OptElement.type == element.mask:
            return np.eye(2)
    if OptElement.type == element.grating:
            return [[1., 0.],[0., 1.]] #DA FARE

def F(yin, z, dndx, dndy, dndz, n): 
    """
    Return derivatives for second-order ODE y'' = - a(z) y' + b(z)
    """
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
    """
    Return derivatives for second-order ODE x'' = - a(z) x' + b(z)
    """
    x=[0, 0]
    x[0]=xin[0] #x(t)
    x[1]=xin[1] #first derivative of x(t)
    
    dxin = [0, 0]    # Create a list to store derivatives.
    dxin[0] = xin[1]    # Store first derivative of x(t).
    dxin[1] = -1*a*x[1]+b    # Store second derivative of x(t).
    return dxin


def Fy(yin, z, a, b): 
    """
    Return derivatives for second-order ODE y'' = - a(z) y' + b(z)
    """
    y=[0, 0]
    y[0]=yin[0] #y(t)
    y[1]=yin[1] #first derivative of y(t)
    
    dyin = [0, 0]    # Create a list to store derivatives.
    dyin[0] = yin[1]    # Store first derivative of y(t).
    dyin[1] = -1*a*y[1]+b    # Store second derivative of y(t).
    return dyin



def RayTracingPropagatorParax(x, y, z, Dx, Dy, Dz, zinit, zend, RayMatrix, n) :
    """
    This function evaluates the ray tracing of the photons of class 'photon' contained in RayMatrix
    propagating in a medium whose refractive index is in matrix n. 
    The rays propagates in the paraxial approximation along the z axis
    
    NB z is the light propagation axis
    
    x, y, z the coordinates of the the matrix n along which the light is propagating
    zinit is the integration start
    zend is the integration end
    RayMatrix is the matrix of ray elements of class 'photon' described in Classes.py
    n is the refractive index matrix  
    THIS VERSION DERIVATES BEFORE PROPAGATING
    """
    import numpy as np
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

def RayTracingPropagatorParax2V(x, y, z, Dx, Dy, Dz, zinit, zend, RayMatrix, n) :
    """
    This function evaluates the ray tracing of the photons of class 'photon' contained in RayMatrix
    propagating in a medium whose refractive index is in matrix n. 
    The rays propagates in the paraxial approximation along the z axis
    
    NB z is the light propagation axis
    
    x, y, z the coordinates of the the matrix n along which the light is propagating
    zinit is the integration start
    zend is the integration end
    RayMatrix is the matrix of ray elements of class 'photon' described in Classes.py
    n is the refractive index matrix  
    THIS VERSION USE THE SPLINE INTERPOLATION TO EVALUATE THE DERIVATIVES
    """
    import numpy as np
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
                x0=[RayMatrix[i][j].x[zcount],  RayMatrix[i][j].dx[zcount]] #initial conditions x and x'
                y0=[RayMatrix[i][j].y[zcount],  RayMatrix[i][j].dy[zcount]] #initial conditions y and y'
                
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
                RayMatrix[i][j].dx[zcount+1]=xsol[1,1]
                RayMatrix[i][j].dy[zcount+1]=ysol[1,1]    
                         
                
    return RayMatrix

def RayTracingPropagatorDUED(x, y, z, RayMatrix, n):
    """
    This function evaluates the ray tracing of the photons of class 'photon' contained in RayMatrix
    propagating in a medium whose refractive index is in matrix n. 
    
    x, y, z the coordinates of the the matrix n along which the light is propagating
    RayMatrix is the matrix of ray elements of class 'photon' described in Classes.py
    n is the refractive index matrix  
    
    More details in DUED plasma code
    """
    import numpy as np
    
    dx=x[1]-x[0]
    dy=y[1]-y[0]
    dz=z[1]-z[0]
    gx, gy, gz=np.gradient(n,[dx,dy,dz])
    

    
    return 0

"""
SCREENS
"""
def printLaserProfile(LasProfile,x,y,z0, savePlot='False', address='none'):
    """
    Print laser profile
    """
    import matplotlib.pyplot as plt
    import pylab as pyl
    
    LasPlot = plt.imshow(LasProfile [:, :, z0], extent=[x[0],x[-1],y[0],y[-1]], cmap=plt.cm.plasma)
    plt.colorbar()
    plt.xlabel('z')
    plt.ylabel('x')
    plt.show
    if savePlot:
        pyl.savefig(address, format='png')
        
    return LasPlot
        
"""
OPTICAL PROPERTIES FUNCTIONS
"""
def Sellmeier(wavelength, medium, ordinary='True', B1=0, B2=0, B3=0, C1=0, C2=0, C3=0):
    """
    This function evaluates the Sellmeier equation for a dispersive medium.
    wavelength in m
    medium is a string with the name of the medium
    for other mediums  RefractiveIndex.info
    LIST OF MEDIA:
        'BK7'
        'fused silica'
        'sapphire', ordinary='True'
        'sapphire', ordinary='False'
    """
    import scipy as sci
    
    if medium=='custom': #it just accept the values in input 
        pass        
    elif medium=='BK7': #SCHOTT glass data sheets
        B1=1.03961212
        B2=0.231792344
        B3=1.01046945
        C1=6.00069867E-3 #um^2
        C2=2.00179144E-2 #um^2
        C3=103.560653 #um^2
    elif medium=='fused silica': #SCHOTT glass data sheets
        B1=0.696166300
        B2=0.407942600
        B3=0.897479400
        C1=4.67914826E-3 #um^2
        C2=1.35120631E-2 #um^2
        C3=97.9340025 #um^2
    elif medium=='sapphire': #SCHOTT glass data sheets
        if ordinary=='True': #or p polarized
            B1=1.43134930
            B2=0.65054713
            B3=5.3414021
            C1=5.2799261E-3 #um^2
            C2=1.42382647E-2 #um^2
            C3=325.017834 #um^2
        elif ordinary=='False': #or s polarized
            B1=1.5039759
            B2=0.55069141
            B3=6.5927379
            C1=5.48041129E-3 #um^2
            C2=1.47994281E-2 #um^2
            C3=402.89514 #um^2
        
    l=wavelength*1E6 #convert the wavelength in 
    n=sci.sqrt(1+(B1*l**2/(l**2-C1))+(B2*l**2/(l**2-C2))+(B3*l**2/(l**2-C3)))
    return n

def ConvertPlasma2n(n_e, LambdaProbe):
    """
    This function converts plasma density in refrective index n
    
    n_e free electron plasma density in [cm-3]
    LambdaProbe photon wavelength in [m]
    n refractive index [adim]
    """
    import scipy as sci
    from scipy.constants import c,e, m_e, epsilon_0, pi
#    c=299792458 #m/s
#    e=1.60219E-19; #C
#    m_e=9.1091E-31; #massa a riposo elettrone	kg 
#    mp=1.6725E-27; #massa a riposo del protone	kg 
#    kb=1.3806488E-23; #costante di Boltzmann J/K
#    eps0=8.854187817e-12; #F/m
    
    n_c=m_e*epsilon_0*((2*pi*c)/e)**2/(LambdaProbe**2)*1e-6 #[cm^-3]
    n=sci.sqrt(1-n_e/n_c)
    return n

"""
LIGHT UTILITIES
"""

def LaserEnergy(LaserProfile, x,y,z):
    """
    It evaluates the energy of the fasor matrix LaserProfile and return its value in [J]
    
    LaserProfile [sqrt(W)/m] (it is the fasor matrix Upb=[Upbxy, Upbz])
    x,y,z in [m]
    """
    from scipy.constants import c
    
    Energy=LaserProfile.sum()/(abs(z[-1]-z[0])/c)*(abs(x[-1]-x[0]))*(abs(y[-1]-y[0]))
    
    return Energy

def Irradiance(U):
    """
    It convert U in irradiance I [W/m^2]
    U fasor of the source field
    """
    import numpy as np
    
    I=np.power(abs(U), 2)
    return I
from __future__ import division
import numpy as np
import math
import utils as ut
from math import pi
from numpy import inf
from scipy import constants as cns
from scipy.special import exp1 as E1
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline
from utils import almost_eq,memoize

Einf=cns.physical_constants ['Rydberg constant times hc in J'][0]
a0=cns.physical_constants ['Bohr radius'][0]*100 #in cm
C0=cns.c/1000 #C0 in km/s
ev=cns.electron_volt
LTE10K=2*np.arange(1,150)**2/np.exp(Einf*(-1.0/np.arange(1,150)**2)/(cns.k*1e4)) # Thermal distribution of electron levels fo 10kK H
LTE10K[49:]=0                                                                    # n=50 has a frequency of 0.1mm assume wavelengths longer than this are optically thin and can be whotoionised by stellar radiation
LTE10K/=LTE10K.sum()                                                             # Normalised
LTE13ev=2*np.arange(1,350)**2/np.exp(Einf*(-1.0/np.arange(1,350)**2)/(Einf))     #Thermal distribution of electron levels fo 20kK H, used for just recombined H atoms, assuming H+ and e- are both thermalised at 10kK
LTE13ev[49:]=0                                                                   # n=50 has a frequency of 0.1mm assume wavelengths longer than this are optically thin and can be whotoionised by stellar radiation
LTE13ev/=LTE13ev.sum()                                                           # Normalised


def boundFreeCross (n, k):
    r"""gives the bound free cross section in units of square Bohr radii for 
the capture to the nth level at normalised wavenumber k 
k=\omega_n / \omega where \omega_n is the frequency of the inf->n transition"""
    if k<=1:return 0.2824407769123843*n*k**3 #64pi/(3^1.5*137)
    else :return 0

nus=np.logspace(8,15,1000)
absorb=np.zeros_like(nus)
for i in xrange(1000):
   for n in xrange(1,150):
        nu_n=Einf*(1.0/n**2)/cns.h
        if nu_n<nus[i]:
            k=nu_n/(nus[i])
            absorb[i]+=(boundFreeCross(n,k)*a0**2) * LTE10K[n-1]

absorb=UnivariateSpline(nus,absorb,k=1,s=0)

gamma=np.array([0,  4.69669480e+08,   9.98102940e+07,   3.01775990e+07,    1.15506010e+07,   5.18981800e+06,   2.90148700e+06,      1.67727900e+06,   1.03433400e+06]) #natural widths of hydrogen levels
_A={}
dat= np.loadtxt('/localhome/pytd/.local/lib/python2.7/site-packages/FreeFree/einsteinA.dat')
k=0
for i in xrange(1,10):
    for j in xrange(i+1,11):
        _A[(j,i)]=dat[k,2]
        k+=1
del dat, k

def einsteinA(u,l):
    assert u>l
    if l>50 and u==l+1:
        return 5.3e9*l**-5
    elif u<=10:
        return _A[(u,l)]
    else:
        raise ValueError("No coefficients for these transitions")
    
def einsteinB(FROM,TO):
    if FROM>TO: u,l=FROM,TO
    else:       l,u=FROM,TO
    B=einsteinA(u,l)*cns.c**2/(2*cns.h*(Einf*(1.0/l**2-1.0/u**2)/cns.h)**3)
    if FROM<TO: B*=u**2/float(l**2)
    return B

def voigt(a,u):
    if np.ndarray in [type(a),type(u)]:
        out=np.zeros_like(a);i=0
        for aa,uu in zip(a.flat,u.flat):
            f=lambda y: np.exp(-y**2)/(aa*aa+(uu-y)**2)
            out.flat[i]= aa/pi * quad(f,-inf, +inf)[0]
            i+=1
        return out
    else:
        f=lambda y: np.exp(-y**2)/(a*a+(u-y)**2)
        return a/pi * quad(f,-inf, +inf)[0]

<<<<<<< HEAD
def _lineAbs_cgs(nu, Ne, Nneut, u, l, T, v=0, lineProfile=1):
    "absorption coeficient due to bound free transitions with the absorbing material moving at v relative to the observer (+ve = blue shifted)"
#    a=0
#    if almost_eq(T, 1e4): 
    N=LTE10K*Nneut
    for uu in xrange(l,150):
        N[uu-1]=+Ne**2*thermalRecomb(uu,T)
    Eline=Einf*(1.0/l**2-1.0/u**2)
    a=absorb([nu*(1+v/C0)])[0]
    if u!=np.inf:
        a+=Eline/(4*pi)*\
          N[l-1]*einsteinB(l,u)*(1-l**2*N[u-1]/(u**2*N[l-1]))
    if lineProfile:
        a*=lineProfile_cgs(nu,u,l,Ne,T,v)
    return a

def lineAbs_cgs(nu, arr, u, l, T, dv=50.0/300000):
    """wrapper for _lineEmmis_cgs
Assumes arr will be a 2x1D vector with arr[0,:]=Ne arr[1,:]=Nneut, arr[2,:]=v"""
#    f=lambda x : _lineAbs_cgs(nu,x[0],x[1],u,l,T,x[2], lineProfile=0)
    f=lambda x : quad(_lineAbs_cgs, nu*(1-dv/2), nu*(1+dv/2), (x[0],x[1],u,l,T,x[2]))[0]
    return np.apply_along_axis(f,0,arr)

def lineProfile_cgs(nu, u,l, Ne, T, v=0):
    nu0=Einf*(1.0/l/l-1.0/u/u)/cns.h
    #nud=dopplerWidth(nu0,T)
    nud=0.0004285850633666147*nu0*math.sqrt(T)
    nus=nu*(1-v/C0)
    if abs(nus-nu0)>5*nud: return 0
    else:                  return math.exp(-(nu*(1-v/C0)-nu0)**2/nud**2) /nud/math.sqrt(pi)   
#ground state has no intrinsic width and unbound has practically none so we run into numerical issues trying to calc a voigt profile, instead just do Doppler
#    if 1:
#        nu_col=4.8e-7*Ne/(T*cns.Boltzmann/ev)**1.5 #ion collision frequency s^-1
#        lam=1.579e5/T
#        a=(gamma[l-1]+ 2 *nu_col)/(4*pi*Vd) 
#        u=(nu-nu0)/Vd
#        return voigt(a,u)/(Vd*math.sqrt(pi))
    

def _thermalRecomb (n,T):
    return 3.262e-6*M(n,T)
thermalRecomb=memoize(_thermalRecomb)

def M(n,T):
    X=Einf*(1.0/n**2)/(cns.Boltzmann*T) #ionisation energy of nth level in units of k_B T
    if X<100:
        return E1(X)*math.exp(X)/(n*math.sqrt(T))**3
    else:
        return (1.0/X - 1.0/X**2 + 2.0/X**3 - 6.0/X**4 + 24.0/X**5) / (n*math.sqrt(T))**3 #n->inf expansion of e^x . E1(x) accurate to ~1+1^-8 at X=100

def _lineEmiss_cgs (nu, Ne, Nn, u, l, T, v=0, lineProfile=1):
    """gives the emissivity at the line centre
ne :elsectron number density (cm^-3)
nu: recombinations to which hydrogen state
T: temperature (K)
v: velocity of the emitting material (+ve= towards observer)"""
    if u==np.inf:   
        nu0=(1+v/C0)*(Einf/l**2)/cns.h
        if abs(nu-nu0)/nu0 < 0.1:
            X=Ne**2*thermalRecomb(l,T)*Einf*(1.0/l**2)*1e7 /(4*pi)
    else:
        X=Einf*(1.0/l**2-1.0/u**2)/(4*pi)*einsteinA(u,l)*\
          (Ne**2*thermalRecomb(u,T)+Nn*LTE10K[u-1])
    if lineProfile :
        X*=lineProfile_cgs(nu,u,l,Ne,T,v)

    return X


def lineEmiss_cgs(nu, arr, u, l, T, dv=50.0/300000):
    """wrapper for _lineEmmis_cgs
Assumes arr will be a 2x1D vector with arr[0,:]=Ne arr[1,:]=v"""
    f=lambda x : quad(_lineEmiss_cgs, nu*(1-dv/2), nu*(1+dv), (x[0],x[1],u,l,T,x[2]))[0]
#    f=lambda x : _lineEmiss_cgs(nu, x[0],x[1],u,l,T,x[2], lineProfile=1)
    return np.apply_along_axis(f,0,arr)

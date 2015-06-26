import numpy as np
import math
from math import pi
from numpy import inf
from scipy import constants as cns
from scipy.special import exp1 as E1
from scipy.integrate import quad
from utils import almost_eq

Einf=cns.physical_constants ['Rydberg constant times hc in J'][0]
C0=cns.speed_of_light
a0=cns.physical_constants ['Bohr radius'][0]*100 #in cm
ev=cns.electron_volt
LTE10K=2*np.arange(1,50)**2/np.exp(-Einf*(1.0/np.arange(1,50)**2)/(cns.Boltzmann*1e4)) # Thermal distribution of electron levels fo 10K H
LTE10K/=LTE10K.sum()                                                                   # Normalised

gamma=np.array([0,  4.69669480e+08,   9.98102940e+07,   3.01775990e+07,    1.15506010e+07,   5.18981800e+06,   2.90148700e+06,      1.67727900e+06,   1.03433400e+06]) #natural widths of hydrogen levels

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

def boundFreeCross (n, k):
    r"gives the bound free cross section in units of square Bohr radii for the capture to the nth level at normalised wavenumber k (k=\omega_n / \omega where \omega_n is the frequency of the inf->n transition"
    if type(k)==np.ndarray:
        return 0.2824407769123843*n*k**3 * (k<=1)
    elif k<=1:
        return 0.2824407769123843*n*k**3 #64pi/(3^1.5*137)
    else :
        return 0

def lineAbs_cgs(nu, Nneut, T, v=0):
    "absorption coeficient due to bound free transitions with the absorbing material moving at v relative to the observer (+ve = blue shifted)"
    a=0
#    if almost_eq(T, 1e4): 
    N=LTE10K
#    else:
#        N=2*np.arange(1,20)**2/np.exp(Einf*(0-1.0/np.arange(1,20)**2)/(cns.Boltzmann*T)) # Thermal distribution of electron levels fo 10K H
    N/=N.sum()                                                                       # Normalised
    for n in xrange(1,10):
        omega_n=Einf*(1.0/n**2)/cns.hbar*(1+v/C0)
        k=omega_n/(2*pi*nu)
        a+=(boundFreeCross(n,k)*a0**2) * N[n-1]
    return a*Nneut

def dopplerWidth(nu0, T):
    return nu0/C0*np.sqrt(2*cns.Boltzmann*T/cns.m_p)

def lineProfile_cgs(nu, n, Ne, T, v=0):

    nu0=(1+v/C0)*(Einf/n**2)/cns.h
    Vd=dopplerWidth(nu0,T)*100 # v in cm/s
#    if n==1:
    return np.exp(-(nu-nu0)**2/Vd**2)/(Vd*math.sqrt(pi)) #ground state has no intrinsic width and unbound has practically none so we run into numerical issues trying to calc a voigt profile, instead just do Doppler
#    else:
#        nu_col=4.8e-7*Ne/(T*cns.Boltzmann/ev)**1.5 #ion collision frequency s^-1
#        lam=1.579e5/T
#        a=(gamma[n-1]+ 2 *nu_col)/(4*pi*Vd) 
#        u=(nu-nu0)/Vd
#        return voigt(a,u)/(Vd*math.sqrt(pi))
    

def thermalRecomb (n,T):
    return 3.262e-6*M(n,T)

def M(n,T):
    X=Einf*(1.0/n**2)/(cns.Boltzmann*T) #ionisation energy of nth level in units of k_B T
    if type(X)==np.ndarray:
        return (E1(X)*np.exp(X)*(X<100)+
                (X>100)*(1.0/X - 1.0/X**2 + 2.0/X**3 - 6.0/X**4 + 24.0/X**5))/(n*np.sqrt(T))**3
    if X<100:
        return E1(X)*math.exp(X)/(n*math.sqrt(T))**3
    else:
        return (1.0/X - 1.0/X**2 + 2.0/X**3 - 6.0/X**4 + 24.0/X**5) / (n*math.sqrt(T))**3 #n->inf expansion of e^x . E1(x) accurate to ~1+1^-8 at X=100

def lineEmiss_cgs (nu, Ne,  T, v=0):
    """gives the emissivity at the frequency nu
ne :elsectron number density (cm^-3)
nu: recombinations to which hydrogen state
T: temperature (K)
v: velocity of the emitting material (+ve= towards observer)"""
    if type(Ne)==np.ndarray:
        X=np.zeros_like(Ne)
        for n in range(1,50):
            nu0=(1+v/C0)*(Einf/n**2)/cns.h
            X+=(lineProfile_cgs(nu,n,Ne,T)*Ne**2*thermalRecomb(n,T)*Einf*(1.0/n**2) * 1e7 /(4*pi))*(abs(nu-nu0)/nu0 < 0.1)
    else:
        X=0
        for n in range(1,50):
            nu0=(1+v/C0)*(Einf/n**2)/cns.h
            if abs(nu-nu0)/nu0 < 0.1:
                X+=lineProfile_cgs(nu,n,Ne,T)*Ne**2*thermalRecomb(n,T)*Einf*(1.0/n**2) * 1e7 /(4*pi)
    return X

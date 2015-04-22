import numpy as np
import math
from math import pi
from scipy import constants as cns
from scipy.special import exp1 as E1
from scipy.integrate import quad

Rinf=cns.physical_constants ['Rydberg constant times hc in J'][0]
C0=cns.speed_of_light
a0=cns.Bohr*100 #in cm
ev=cns.electron_volt

def voigt(a,u):
    f=lambda y: exp(-y**2)/(a*a+(u-y)**2)
    return a/pi * quad(f,-inf, +inf)[0]


def boundFreeCross (n, k):
    r"gives the bound free cross section in units of square Bohr radii for the capture to the nth level at normalised wavenumber k (k=\omega_n / \omega where \omega_n is the frequency of the inf->n transition"
    if k<=1:return 64*np.pi/(3**1.5*137)*n*k**3
    else :return 0

def lineAbs_cgs(nu, Nn, v=0):
    "absorption coeficient due to bound free transitions with the absorbing material moving at v relative to the observer (+ve = blue shifted)"
    a=0
    for n in xrange(20):
        omega_n=Rinf*(1.0/n**2)/cns.hbar*(1+v/C0)
        k=omega_n/(2*pi*nu)
        a+=boundFreeCross(n,k)*a0**2
    return a*Nn

def dopplerWidth(nu0, T):
    return nu0/C0*math.sqrt(cns.Boltzmann*T/cns.m_p)

def lineProfile_cgs(nu, nu0, Ne, T):
    nu_col=4.8e-7*Ne/(T*cns.Boltzmann/ev)**1.5 #ion collision frequency s^-1
    Vd=dopplerWidth(nu0,T)
    a=2*nu_col/(4*pi*Vd)
    u=(nu-nu0)/Vd
    return voigt(a,u)/(Vd*sqrt(pi))
    
def thermalRecomb (n,T):
    return 3.262e-6*M(n,T)

def M(n,T):
    X=Rinf*(1.0/n**2)/(cns.Boltzmann*T) #ionisation energy of nth level in units of k_B T
    
    if X<100:
        return E1(X)*math.exp(X)/(n*math.sqrt(T))**3
    else:
        return (1.0/X - 1.0/X**2 + 2.0/x**3 - 6.0/x**4 + 24.0/x**5) / (n*math.sqrt(T))**3 #n->inf expansion if e^x . E1(x) accurate to ~1+1^-8 at X=100

def lineEmiss_cgs (Ne, n, T):
    """gives the emissivity at the line centre
ne :elsectron number density (cm^-3)
n: recombinations to which hydrogen state
T: temperature (K)"""
    return Ne**2*thermalRecomb(n,T)*Rinf*(1.0/n**2) * 1e7 /(4*pi)


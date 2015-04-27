from numpy import *
from gaunt import gaunt
from scipy import constants as cns

def emiss(npls,t , nu):
    "free free emisivity, in CGS unit (blegh)"
    
    g=[gaunt(t,z,nu) for z in [1,2,3,6]]# gaunt factors are of order unity 
    npl1,npl2,npl3,npl6=npls

    ne= 2*npl2 
    ne+= npl1
    ne+= 3*npl3 
    ne+= 6*npl6

    ni=npls.sum(0)
    hnukt = cns.Planck*nu/(cns.Boltzmann*t) #(h nu) / (k_b T)
    gfac=(g[0]*npl1+g[1]*4.0*npl2+g[2]*9.0*npl3+g[3]*36.0*npl6) # g_nu*Z**2
    ni*=gfac

    epsff = exp(-hnukt) # ff emission co-efficient from rybiki & lightman 5.14b
    epsff*=ni
    epsff/=sqrt(t)
    epsff*=6.8e-38/(4*pi)

    kapff = 0.01765/nu # ff absorption co-efficient from rybiki & lightman 5.19b
    kapff/=nu
    kapff/=pow(t,1.5)
    kapff*=ni
    
    return epsff,kapff,g

def eDensity(rho,t):

    nh  =rho/(cns.m_p*1000)*0.89
    nhe =rho/(cns.m_p*1000*4)*0.1
    ncno=rho/(cns.m_p*1000*14.24)*0.01

    npl1=nh+(t<3e5)*nhe
    npl2=(t<3e5)*ncno+(t>3e5)*nhe
    npl3=logical_and(t>3e5,t<1e6)*ncno
    npl6=(t>1e6)*ncno

    return  array([npl1, npl2, npl3,  npl6])

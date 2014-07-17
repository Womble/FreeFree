from numpy import *
from gaunt import gaunt
from scipy import constants as cns

def emiss(npls,t , nu):
    "free free emisivity, in CGS unit (blegh)"
    
    g=[gaunt(t,z,nu) for z in [1,2,3,6]]# gaunt factors are of order unity 
    npl1,npl2,npl3,npl6=npls
    ne=npl1 + 2*npl2 + 3*npl3 + 6*npl6
    hnukt = 4.801449e-11*nu/t #(h nu) / (k_b T)
    epsff = (6.8e-38*exp(-hnukt)*ne*(g[0]*npl1+g[1]*4.0*npl2+g[2]*9.0*npl3+g[3]*36.0*npl6))/sqrt(t) # ff emission co-efficient from rybiki & lightman 5.14b
    epsff /=(4*pi)
    kapff = (0.01765/((nu*nu)*pow(t,1.5)))*ne*(g[0]*npl1 + g[1]*4.0*npl2 + g[2]*9.0*npl3) # ff absorption co-efficient from rybiki & lightman 5.19b

    return epsff,kapff

def eDensity(rho,t):

    nh  =rho/(cns.m_p*1000)*0.89
    nhe =rho/(cns.m_p*1000*4)*0.1
    ncno=rho/(cns.m_p*1000*14.24)*0.01

    npl1=nh+(t<1e6)*nhe
    npl2=(t<3e5)*ncno+(t>3e5)*nhe
    npl3=logical_and(t>3e5,t<1e6)*ncno
    npl6=(t>1e6)*ncno

    return  [npl1, 2.0*npl2, 3.0*npl3,  6.0*npl6]

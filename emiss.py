from numpy import *
from gaunt import gaunt
from scipy import constants as cns
from scipy.interpolate import UnivariateSpline
from utils import powerfit

_dOs=UnivariateSpline(*zip(*[
    [1.000E+00, 10410.0],
    [1.166E+00 ,8915.0],
    [1.359E+00 ,7704.0],
    [1.585E+00 ,6541.0],
    [1.848E+00 ,5595.0],
    [2.154E+00 ,4756.0],
    [2.512E+00 ,3968.0],
    [2.712E+00 ,3598.0],
    [2.818E+00 ,7410.0],
    [2.929E+00 ,31380.0],
    [3.043E+00 ,35060.0],
    [3.162E+00 ,19430.0],
    [3.286E+00 ,6221.0],
    [3.414E+00 ,3694.0],
    [3.548E+00 ,3067.0],
    [3.687E+00 ,2224.0],
    [3.831E+00 ,2169.0],
    [3.981E+00 ,2043.0],
    [4.299E+00 ,2073.0],
    [4.642E+00 ,2205.0],
    [5.012E+00 ,1497.0],
    [5.412E+00 ,1353.0],
    [5.843E+00 ,2876.0],
    [6.310E+00 ,2579.0],
    [6.813E+00 ,2386.0],
    [7.356E+00 ,1548.0],
    [7.943E+00 ,1206.0],
    [8.577E+00 ,1107.0],
    [8.912E+00 ,1744.0],
    [9.261E+00 ,1696.0],
    [9.623E+00 ,2773.0],
    [9.810E+00 ,5478.0],
    [1.000E+01 ,3138.0],
    [1.039E+01 ,3479.0],
    [1.080E+01 ,5101.0],
    [1.122E+01 ,7062.0],
    [1.166E+01 ,8523.0],
    [1.259E+01 ,8866.0],
    [1.359E+01 ,7574.0],
    [1.468E+01 ,5813.0],
    [1.585E+01 ,4628.0],
    [1.711E+01 ,4059.0],
    [1.848E+01 ,3698.0],
    [1.995E+01 ,2996.0],
    [2.155E+01 ,2271.0],
    [2.512E+01 ,1248.0],
    [2.929E+01 ,584.0],
    [3.415E+01 ,835.6],
    [3.687E+01 ,1084.0],
    [3.981E+01 ,1362.0],
    [4.299E+01 ,1633.0],
    [4.642E+01 ,1653.0],
    [5.412E+01 ,1191.0],
    [6.310E+01 ,820.7],
    [7.357E+01 ,540.9],
    [8.577E+01 ,339.6],
    [1.000E+02 ,171.5],
    [1.166E+02 ,113.4],
    [1.359E+02 ,79.15],
    [1.585E+02 ,55.72],
    [1.848E+02 ,38.42],
    [2.260E+02 ,25.09],
    [3.500E+02 ,11.22],
    [5.000E+02 ,5.502],
    [7.000E+02 ,2.815],
    [1.000E+03 ,1.481],
    [1.300E+03 ,0.9621]]), k=1, s=0)
_dOcoeff=powerfit(logspace(2.7,log10(1300),500),_dOs(logspace(2.7,log10(1300),500)))

def dustOpacity(Lambda):
    """returns linearly interpolated dust opacities (in cm^-2 / g)
from wavelengths (in microns) between 1 and 1300 micron.
For wavelengths greater than 1300 micron uses a power law fited to the data beween 500 and 1300 micron"""
    assert Lambda>1
    if Lambda<=1300:
        return _dOs(Lambda)
    else:
        return _dOcoeff[1]*Lambda**_dOcoeff[0]
    

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

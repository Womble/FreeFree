from numpy import *
from scipy import constants as cns
from emiss import emiss,eDensity
from scipy.ndimage.interpolation import rotate

from ctypes import *
from numpy.ctypeslib import ndpointer
try:
    libpath='/home/student13/pytd/home.linux/.local/lib/python2.7/site-packages/FreeFree/integ.so'
    lib=CDLL(libpath)
    func=lib.integrate 
    print 'importing integration routine from %s'%libpath
    func.argtypes = [ndpointer(c_double), c_double, c_size_t]
    func.restype = c_double
except:
    print "cant import c funtion for integrating from "+libpath+", has integ.c been compiled with -fPIC and -shared?"

PC2CM = 3.085678e18 #1pc in cm
SQRAD2STR = 1/(pi*0.25) #convert square radians to steradians

def integratePY (column, cellLength):
    "columns should be a 1D array containing dt, epsilons and kappas  in the first second and third third /// this is the heavy lifting part and should be ported to C"
    intensity=0
    l=column.size
    if l%2!=0 : raise ValueError('column cannont be split into halves')
    for i in xrange(l/2):
        n=max(1,int(2*column[2*i]))
        dt=column[2*i]/n
        for j in xrange(n):
            intensity+=dt*(column[2*i+1]-intensity)
    return intensity

def doppler (vr):
    "dopller shift, +ve= towards observer, assumes given in lightspeed units if no values >1 else in cm/s"
    if (abs(vr)>=1).any(): vr/cns.speed_of_light/100
    return sqrt((1+vr)/(1-vr))

class freeFree():
    def __init__(self,Rho, Temp, Length, v=None):
        "Rho in g/cm^3"
        self.rho=Rho
        self.t=Temp
        self.npls=eDensity(Rho,Temp)
        self.length=Length #length per unit cell (in cms, ew)

    def ne(self):
        return self.npls[0]+self.npls[1]*2+self.npls[2]*3+self.npls[3]*6
    
    def epsNkap(self,nu):
        self.eps,self.kap,self.gaunts=emiss(self.npls, self.t, nu)
        #assume emission lines are optically thin themselves so no absorbtion coefficient from photo-ionoisation
        #centre_nu=nu0*doppler(self.v[...,2])
        #sigma2=(cns.Boltzmann*self.t/cns.m_p/cns.speed_of_light**2)*centre_nu**2
        #line_eps+=self.ne()*self.npls[0]* 2.076e-11*2.2/np.sqrt(self.t) * exp(-(nu-centre_nu)**2/2/sigma2)

    def taus(self,nu):
        self.epsNkap(nu)
        self.dt=self.kap*self.length

    def rotatecube(self,theta=0,phi=0):
        rho=self.rho
        temp=self.t
        if (int(phi)%360)!=0:
            rho =rotate(self.rho,phi,  (0,1), mode='nearest', order=1)
            temp=rotate(self.t,phi, (0,1), mode='nearest', order=1)
        if (int(theta)%360)!=0:
            rho =rotate(rho,theta, (1,2), mode='nearest', order=1)
            temp=rotate(temp,theta,(1,2), mode='nearest', order=1)
        rho[rho<0]=1e-30
        temp[temp<1]=1
        self.rho=rho
        self.t=temp
        
    def rayTrace(self,nu,theta=0,phi=0, dist=500, returnRotatedCube=0):
        "integrate along the specified axis after rotating the cube through phi and theta (in deg)"
        if dist<1e9: dist*=PC2CM #assume distances less than 10^9 are given im parsecs, larger in cm

        try:
            flag=self.dt.any()
        except:
            flag=0
        if theta==0 and phi==0 and flag:
            tempcube=self
            print 'reusing dt'

        else :
            rho=self.rho
            temp=self.t
            if (int(phi)%360)!=0:
                rho =rotate(rho,phi,  (0,1), mode='nearest', order=1)
                temp=rotate(temp,phi, (0,1), mode='nearest', order=1)
                #need to include rotation to vel field in here              
            if (int(theta)%360)!=0:
                rho =rotate(rho,theta, (1,2), mode='nearest', order=1)
                temp=rotate(temp,theta,(1,2), mode='nearest', order=1)
            rho[rho<0]=1e-30
            temp[temp<1]=1
            tempcube=freeFree(rho, temp, self.length)
            tempcube.taus(nu)
        f=lambda x : func(x,self.length, x.size)
#        f=lambda x : integratePY(x,self.length)
        s=tempcube.dt.shape
        source=(tempcube.eps/tempcube.kap)
        source[isnan(source)]=0
        arr=empty((s[0],s[1],2*s[2]))
        tempcube.dt[tempcube.dt<1e-30]=1e-30 #dont allow tau of cell to be less than 1e-30
        arr[:,:,::2]=tempcube.dt[::-1]   #invert z axis so we integrate from the far side towards observer over the z axis
        arr[:,:,1::2]=source[::-1]       #source function (eps/kap)
        out=apply_along_axis(f,-1,arr)    # to c after the rotation and do the integrating there
        #apply_along_axis takes 1D z slices through the cube and passes them to integrate
        pix =abs(self.length/dist)
        print pix
        if returnRotatedCube:return out*pix*pix*SQRAD2STR*1e26,tempcube 
        else :               return out*pix*pix*SQRAD2STR*1e26 #output in mJy/pix

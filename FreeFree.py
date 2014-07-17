from numpy import *
from scipy import constants as cns
from emiss import emiss,eDensity
from scipy.ndimage.interpolation import rotate

from ctypes import *
from numpy.ctypeslib import ndpointer
try:
    libpath='/localhome/pytd/opt/radio-py/integ.so'
    lib=CDLL(libpath)
    func=lib.integrate # for example
    func.argtypes = [ndpointer(c_double), c_double, c_size_t]
    func.restype = c_double
except:
    print "cant import c funtion for integrating from "+libpath+", has integ.c been compiled with -fPIC and -shared?"

PC2CM = 3.085678e18 #1pc in cm
SQRAD2STR = 1/(pi*0.25) #convert square radians to steradians

#def integrate (column, cellLength):
#    "columns should be a 1D array containing dt, epsilons and kappas  in the first second and third third /// this is the heavy lifting part and should be ported to C"
#    intensity=0
#    l=column.size
 #   if l%3!=0 : raise ValueError('column cannont be split into thirds')
 #   for i in xrange(l/3):
 #       expDt=exp(-column[3*i]/cellLength)        
 #       intensity= intensity*expDt+column[3*i+1]*(1-expDt)/column[3*i+2]  
 #   return intensity

class freeFree():
    def __init__(self,Rho, Temp, Length):
        self.rho=Rho
        self.t=Temp
        self.npls=eDensity(Rho,Temp)
        self.length=Length #length per unit cell (in cms, ew)

    def ne(self):
        return self.npls[0]+self.npls[1]*2+self.npls[2]*3+self.npls[3]*6

    def epsNkap(self,nu):
        self.eps,self.kap=emiss(self.npls, self.t, nu)
    
    def taus(self,nu):
        self.epsNkap(nu)
        self.dt=self.kap*self.length
        
    def rayTrace(self,nu,theta=0,phi=0, dist=500):
        "integrate along the specified axis after rotating the cube through phi and theta (in deg)"
        if dist<1.5e13: dist*=PC2CM #assume distances less than 1au are in parsecs, otherwise in cm

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
            if (int(theta)%360)!=0:
                rho =rotate(rho,theta, (1,2), mode='nearest', order=1)
                temp=rotate(temp,theta,(1,2), mode='nearest', order=1)
            tempcube=freeFree(rho, temp, self.length)
            tempcube.taus(nu)
        f=lambda x : func(x,self.length, x.size)
        s=tempcube.dt.shape
        arr=empty((s[0],s[1],2*s[2]))
        arr[:,:,1::2]=tempcube.eps[::-1] #emission,coef * celllength 
        arr[:,:,::2]=tempcube.kap[::-1]   #invert z axis so we integrate from the far side towards observer over the z axis
        out=apply_along_axis(f,-1,arr)    # to c after the rotation and do the integrating there
        #apply_along_axis takes 1D z slices through the cube and passes them to integrate
        pix =abs(self.length/dist)
        return out*pix*pix*SQRAD2STR*1e26 #output in mJy/pix

from numpy import *
from scipy import constants as cns
from scipy.ndimage.interpolation import rotate
from scipy.stats.mstats import gmean
from utils import almost_eq

from emiss import emiss,eDensity


#from ctypes import *
#from numpy.ctypeslib import ndpointer
#try:
#    libpath='/home/student13/pytd/home.linux/.local/lib/python2.7/site-packages/FreeFree/integ.so'
#    lib=CDLL(libpath)
#    Cfunc=lib.integrate 
#    print 'importing integration routine from %s'%libpath
#    Cfunc.argtypes = [ndpointer(c_double), ndpointer(c_size_t), ndpointer(c_double, flags='W')]
#    Cfunc.restype = c_double
#except:
#    print "cant import c funtion for integrating from "+libpath+", has integ.c been compiled with -fPIC and -shared?"

PC2CM = 3.085678e18 #1pc in cm
SQRAD2STR = 4/pi #convert square radians to steradians

def integrate (source, dt):
    assert source.shape == dt.shape
    s=source.shape
    out=zeros_like(source[...,0])
    tmp=out.copy()
    mask=zeros_like(out, dtype=bool)
    for i in xrange(s[-1]):
        if dt[...,i].max()>0.1:
            n=5+min(20, int(0.75/dt[...,i].max()))
            for _ in xrange(n):
                tmp[...]=source[...,i]
                tmp-=out
                tmp*=dt[...,i]
                tmp/=n
                out+=tmp
            mask[...]=dt[...,i]>10
            out[mask]=source[...,i][mask]
        else:
            tmp[...]=source[...,i]
            tmp-=out
            tmp*=dt[...,i]
            out+=tmp
    return out


def doppler (vr):
    "dopller shift, +ve= towards observer, assumes given in lightspeed units if no values >1 else in cm/s"
    if (abs(vr)>=1).any(): vr/cns.speed_of_light/100
    return sqrt((1+vr)/(1-vr))

def trimCube(cube, thresh):
    "returns the silce which trims planes off cube if all the values in the plane are < thresh"
    sl=[]
    for j,s in enumerate(cube.shape):
        a,b,i,flag=0,0,0,0
        while i<s and not(flag):
            if j==0:
                if (cube[i,:,:]<thresh).all(): a+=1
                else                 : flag=1
            elif j==1:
                if (cube[:,i,:]<thresh).all(): a+=1
                else                 : flag=1
            else :
                if (cube[:,:,i]<thresh).all(): a+=1
                else                 : flag=1
            i+=1
        i,flag=0,0
        while i<s and not(flag):
            if j==0:
                if (cube[-i,:,:]<thresh).all(): b-=1
                else                 : flag=1
            elif j==1:
                if (cube[:,-i,:]<thresh).all(): b-=1
                else                 : flag=1
            else :
                if (cube[:,:,-i]<thresh).all(): b-=1
                else                  : flag=1
            i+=1
        sl.append(slice(a,b))
    return sl

class freeFree():
    def __init__(self,Rho, Temp, Length):#, v=None):
        """Rho is ion density cube in g/cm^3
Temp is temperature cube in K (Rho and Temp need to have the same shape)
Length is the size of one cell in the Rho and Temp cubes"""
        self.rho=Rho
        self.t=Temp
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
        self.npls=eDensity(self.rho,self.t)
        self.epsNkap(nu)
        self.dt=self.kap*self.length

    def rotatecube(self,theta=0,phi=0, trim=0):
        "angles in degrees"
        rho=self.rho.copy()
        t=self.t.copy()
        if (int(theta)%360)!=0:
            rho =rotate(rho,theta, (1,2), mode='nearest', order=1)
            t=   rotate(t,  theta, (1,2), mode='nearest', order=1)
        if (int(phi)%360)!=0:
            rho =rotate(rho,phi, (0,1), mode='nearest', order=1)
            t=rotate   (t,  phi, (0,1), mode='nearest', order=1)
        rho[rho<1e-30]=1e-30
        t[t<1]=1
        thresh=rho
        if trim:
            for _ in rho.shape:
                thresh=gmean(thresh)
            sl=trimCube(rho, thresh*5)
            self.rho=rho[sl]
            self.t=t[sl]
        else:
            self.rho=rho
            self.t=t
        try:
            self.dt[...]=0
        except AttributeError:
            None
        
    def rayTrace(self,nu,theta=0,phi=0, dist=500, returnRotatedCube=0, transpose=0):
        "integrate along the specified axis after rotating the cube through phi and theta (in deg)"
        if dist<1e9: dist*=PC2CM #assume distances less than 10^9 are given im parsecs, larger in cm

        try:
            flag=self.dt.any() and almost_eq(nu,self.lastnu)
        except:
            flag=0
        if theta==0 and phi==0:
            tempcube=self
            if flag:
                print 'reusing dt'
            else:
                print 'calculating taus'
                self.taus(nu)
                self.lastnu=nu
        else :
            tempcube=freeFree(self.rho.copy(), self.t.copy(), self.length)
            tempcube.rotatecube(theta,phi)
            print 'calculating taus'
            tempcube.taus(nu)
#        f=lambda x : Cfunc(x,self.length, x.size)
#        f=lambda x : integratePY(x,self.length)
        s=tempcube.dt.shape
        source=(tempcube.eps/(4*pi)/tempcube.kap)
        source[isnan(source)]=0
        tempcube.dt[tempcube.dt<1e-30]=1e-30 #dont allow tau of cell to be less than 1e-30
        print('integrating')
        if transpose:out=integrate((source.T)[...,::-1],(tempcube.dt.T)[...,::-1]) #integrate from back to front so we are looking down from from +z
        else:        out=integrate(source[...,::-1],tempcube.dt[...,::-1]) 
        pix =abs(self.length/dist)
        self.im=out*pix*pix*SQRAD2STR*1e26
        if returnRotatedCube:return self.im,tempcube 
        else :               return self.im #output in mJy/pix


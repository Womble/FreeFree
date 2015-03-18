from numpy import *
from scipy import constants as cns
from scipy.ndimage.interpolation import rotate
from scipy.stats.mstats import gmean
from utils import almost_eq

from emiss import emiss,eDensity

plankCGS=cns.h*1e7
cCGS=cns.c*100
kbCGS=cns.Boltzmann*1e7

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
        source=(tempcube.eps/tempcube.kap)
        source[isnan(source)]=0
        tempcube.dt[tempcube.dt<1e-30]=1e-30 #dont allow tau of cell to be less than 1e-30
        print('integrating')
        if transpose:out=integrate((source.T)[...,::-1],(tempcube.dt.T)[...,::-1]) #integrate from back to front so we are looking down from from +z
        else:        out=integrate(source[...,::-1],tempcube.dt[...,::-1]) 
        pix =abs(self.length/dist)
        self.im=out*pix*pix*SQRAD2STR*1e23*1000
        if returnRotatedCube:return self.im,tempcube 
        else :               return self.im #output in mJy/pix

def bb(T, nu):
    return 2*plankCGS*nu**3/cCGS**2 * 1/(exp(plankCGS*nu/kbCGS/T)-1)


def test():
    "check RT is working by comparing the fluxes from an optically thick and thin sphere to the analytic formulae"
    x,y,z=mgrid[-1:1:100j,-1:1:100j,-1:1:100j]
    rho=ones((100,100,100), dtype=float)*cns.m_p*1000*10
    rho[sqrt(x*x+y*y+z*z)>0.9]*=1.0e-10
    temp=ones_like(rho)*1.0e4
    RT=freeFree(rho,temp,1.5e11) #RT for a sphere 1au in diameter @T=10,000 n=10/cc
    thinim=RT.rayTrace(1.0e9)
    ne=rho/(cns.m_p*1000)
    thinAnalytic=(6.8e-38*ne**2*RT.gaunts[0]/sqrt(1.0e4)*exp(-cns.h*1e9/cns.Boltzmann/1.0e4)).sum()*RT.length**3*1e23/(4*pi*(500*PC2CM)**2)
    assert almost_eq(thinim.sum()/1000,thinAnalytic, diff=0.1)
    print "optically thin test ok (ratio %.3f)"%(thinim.sum()/1000/thinAnalytic) #usually a bit out as analytic only inculdes hydrogen, within 10% at low temps
    rho*=1e7                    #optically thick sphere
    RT=freeFree(rho,temp,1.5e11)
    thickim=RT.rayTrace(1.0e9)
    thickAnalytic=pi*(RT.length*50)**2*bb(1.0e4,1.0e9)*1e23/(500*PC2CM)**2
    assert almost_eq(thickim.sum()/1000,thickAnalytic, diff=0.1)
    print "optically thin test ok (ratio %.3f)"%(thickim.sum()/1000/thickAnalytic)
    return RT

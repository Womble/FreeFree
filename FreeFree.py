from numpy import *
from scipy import constants as cns
from scipy.ndimage.interpolation import rotate
from scipy.stats.mstats import gmean
from utils import almost_eq,Rx,Ry,Rz,Planck
import gc

from emiss import emiss,eDensity, dustOpacity
import line

plankCGS=cns.h*1e7
cCGS=cns.c*100
kbCGS=cns.Boltzmann*1e7

Gas2Dust=100 #neutral gas to dust mass ratio
DustDestTemp=3000 #dust destruction temperature

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
    for i in xrange(s[-1]):
        out=source[...,i]+exp(-dt[...,i])*(out-source[...,i])
    return out


def doppler (vr):
    "dopller shift, +ve= towards observer, assumes given in lightspeed units if no values >1 else in cm/s"
    if (abs(vr)>=1).any(): vr/c0
    return sqrt((1+vr)/(1-vr))

def trimCube(cube, thresh):
    "returns the silce which trims planes off cube if all the values in the plane are < thresh"
    sl=[]
    for j,s in enumerate(cube.shape):
        a,b,i,flag=0,s,0,False
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
    def __init__(self,RhoI, RhoN, Temp,  Velocity, Length):#, v=None):
        """Rho is ion density cube in g/cm^3
Temp is temperature cube in K (Rho and Temp need to have the same shape)
Length is the size of one cell in the Rho and Temp cubes"""
        self.rho=RhoI
        self.rhoN=RhoN
        self.t=Temp
        self.V=Velocity
        self.length=Length #length per unit cell (in cms, ew)
        self.eps=zeros_like(self.rho)
        self.kap=zeros_like(self.rho)

    def ne(self):
        return self.npls[0]+self.npls[1]*2+self.npls[2]*3+self.npls[3]*6
    
    def epsNkap(self,nu):
        self.eps,self.kap,self.gaunts=emiss(self.npls, self.t, nu)

    def cleartau(self):
        del self.dt
        del self.eps
        del self.kap
        
    def taus(self,nu, ff=1, lines=0, dust=0):
        self.npls=eDensity(self.rho,self.t)
        self.eps=zeros_like(self.rho, dtype=float64)
        self.kap=zeros_like(self.rho, dtype=float64)
        if ff: self.epsNkap(nu)
        if lines :
            u,l=lines
            print u,l,nu
            ne= 2*self.npls[1,...] 
            ne+= self.npls[0,...]
            ne+= 3*self.npls[2,...]
            ne+= 6*self.npls[3,...]
            self.kap.flat+=line.lineAbs_cgs(nu, array([ne.flatten(),
                                                       self.rhoN.flatten()/cns.m_p,
                                                       self.V[2,...].flatten()]),u,l, self.t.mean()) #velocity is vz (ie los) in cm/s
            self.eps.flat+=line.lineEmiss_cgs(nu,array([ne.flatten(), 
                                                        self.rhoN.flatten()/cns.m_p, 
                                                        self.V[2,...].flatten()]),u,l,self.t.mean())
        if dust:
            mask=self.t<DustDestTemp
            d_kap=dustOpacity(cns.speed_of_light/nu*1e6)*self.rhoN[mask]/Gas2Dust
            self.kap[mask]+=d_kap
            self.eps[mask]+=d_kap*Planck(nu,self.t[mask],cgs=True)
#        self.kap[self.kap<1e-40]=1e-40
#        self.eps[self.eps<1e-40]=1e-40
        self.dt=self.kap*self.length

    def rotatecube(self,theta=0,phi=0, trim=0):
        "angles in degrees"
        if (int(theta)%360)!=0:
            rho=rotate(self.rho,theta, (0,2), mode='nearest', order=1)
            rhoN=rotate(self.rhoN,theta, (0,2), mode='nearest', order=1)
            t =rotate(self.t,  theta, (0,2), mode='nearest', order=1)
            v =rotate(self.V,  theta, (1,3), mode='nearest', order=1)
            M=Ry(theta*pi/180)
            M[abs(M)<(finfo(1.0).eps*10)]=0 #set numbers with abs value less than 10 times the floating point epsilon to 0
            f=lambda x : (x*M).flat
            v=apply_along_axis(f,0, v)
            if (int(phi)%360)!=0:
                rho=rotate(rho,phi, (0,1), mode='nearest', order=1)
                rhoN=rotate(rhoN,phi, (0,1), mode='nearest', order=1)
                t  =rotate(t,  phi, (0,1), mode='nearest', order=1)
                v  =rotate(v,  phi, (1,2), mode='nearest', order=1)
                M=Rz(phi*pi/180)
                M[abs(M)<(finfo(1.0).eps*10)]=0
                f=lambda x : (x*M).flat
                v=apply_along_axis(f,0, v)
        elif (int(phi)%360)!=0:
            rho=rotate(self.rho,phi, (0,1), mode='nearest', order=1)
            rhoN=rotate(self.rhoN,phi, (0,1), mode='nearest', order=1)
            t  =rotate(self.t,  phi, (0,1), mode='nearest', order=1)
            v  =rotate(self.V,  phi, (1,2), mode='nearest', order=1)
            M=Rz(phi*pi/180)
            M[abs(M)<(finfo(1.0).eps*10)]=0
            f=lambda x : (x*M).flat
            v=apply_along_axis(f,0, v)
        else:
            rho=self.rho.copy()
            t=self.t.copy()
            v=self.V.copy()
        rho[rho<1e-30]=1e-30
        t[t<1]=1
        thresh=rho
        if trim:
            for _ in rho.shape:
                thresh=gmean(thresh, axis=0)
            sl=trimCube(rho, thresh*5)
            self.rho=rho[sl]
            self.rhoN=rhoN[sl]
            self.t=t[sl]
            sl=[slice(0,3)]+sl
            self.V=v[sl]
        else:
            self.rho=rho
            self.rhoN=rhoN
            self.t=t
            self.V=v
        try:
            self.dt[...]=0
        except AttributeError:
            None
        
    def rayTrace(self,nu,theta=0,phi=0, dist=500, ff=1, lines=0, returnRotatedCube=0, transpose=0, suppressOutput=False):
        "integrate along the specified axis after rotating the cube through phi and theta (in deg)"
        if dist<1e9: dist*=PC2CM #assume distances less than 10^9 are given im parsecs, larger in cm

        try:
            flag=self.dt.any() and almost_eq(nu,self.lastnu)
        except:
            flag=0
        if theta==0 and phi==0:
            tempcube=self
            if flag:
                if not(suppressOutput): print 'reusing dt'
            else:
                if not(suppressOutput): print 'calculating taus'
                self.taus(nu,ff,lines)
                self.lastnu=nu
        else :
            tempcube=freeFree(self.rho.copy(), self.rhoN.copy(), self.t.copy(), self.V.copy(), self.length)
            tempcube.rotatecube(theta,phi)
            if not(suppressOutput): print 'calculating taus'
            tempcube.taus(nu)
#        f=lambda x : Cfunc(x,self.length, x.size)
#        f=lambda x : integratePY(x,self.length)
        s=tempcube.dt.shape
        source=(tempcube.eps/tempcube.kap)
        source[source!=source]=0
        tempcube.dt[tempcube.dt<1e-30]=1e-30 #dont allow tau of cell to be less than 1e-30
        if not(suppressOutput): print('integrating')
        if transpose:out=integrate((source.T)[...,::-1],(tempcube.dt.T)[...,::-1]) #integrate from back to front so we are looking down from from +z
        else:        out=integrate(source[...,::-1],tempcube.dt[...,::-1]) 
        pix =abs(self.length/dist)
        self.im=out*pix*pix*SQRAD2STR*1e23*1000
        if returnRotatedCube:return self.im,tempcube 
        else :               return self.im #output in mJy/pix

    def spectrum(self,nus, theta=0, phi=0, dist=500, trim=1):
        if theta or phi:
            self.rotateCube(theta, phi, trim)
        vals=[]
        for nu in nus:
            try: 
                del self.dt
                del self.npls
                del self.eps
                del self.kap
                del self.gaunts
            except AttributeError:
                None
            gc.collect()
            if nu<1e6 : nu*=1e9 # assume vals < 1MHz are intended to be in GHz
            vals.append(self.rayTrace(nu, dist=dist, suppressOutput=True).sum())
        return vals

def bb(T, nu):
    return 2*plankCGS*nu**3/cCGS**2 * 1/(exp(plankCGS*nu/kbCGS/T)-1)


def test():
    "check RT is working by comparing the fluxes from an optically thick and thin sphere to the analytic formulae"
    x,y,z=mgrid[-1:1:100j,-1:1:100j,-1:1:100j]
    rho=ones((100,100,100), dtype=float)*cns.m_p*1000*1
    rho[sqrt(x*x+y*y+z*z)>0.999]*=1.0e-10
    temp=ones_like(rho)*1.0e4

    RT=freeFree(rho,zeros_like(rho),temp, zeros_like([rho,rho,rho]),1.5e11) #RT for a sphere 1au in diameter @T=10,000 n=1/cc
    thinim=RT.rayTrace(2.0e9)
    ne=rho/(cns.m_p*1000)
    thinAnalytic=(6.8e-38*ne**2*RT.gaunts[0]/sqrt(1.0e4)*exp(-cns.h*2e9/cns.Boltzmann/1.0e4)).sum()*RT.length**3*1e23/(4*pi*(500*PC2CM)**2)
    assert almost_eq(thinim.sum()/1000,thinAnalytic, diff=0.1)
    print "optically thin test ok (ratio %.3f)"%(thinim.sum()/1000/thinAnalytic) #usually a bit out as analytic only inculdes hydrogen, within 10% at low temps

    thinPow=[log10(RT.rayTrace(x*1e8, suppressOutput=True).sum()) for x in xrange(1,11)]
    thinpf=polyfit([log10(x) for x in xrange(1,11)], thinPow,1)[0]
#    assert almost_eq(thinpf, -0.1, 0.05)
    print "powerlaw spectrum for thin sphere",thinpf

    #thin line test
    RT.cleartau()
    Bgam=138475490718753.05
    lineim=RT.rayTrace(Bgam, ff=0,lines=(7,4))
    lineAnalytic=line.Einf*(4.0**-2-7.0**-2)*line.LTE10K[7-1]*\
                 RT.rho.sum()*RT.length**3/cns.m_p * line.einsteinA(7,4)*\
                 1e23/(4*pi*(500*PC2CM)**2)

    print lineim.sum(), lineAnalytic

    rho*=1e8                    #optically thick sphere
    RT=freeFree(rho,zeros_like(rho),temp,zeros_like([rho,rho,rho]),1.5e11)
    thickim=RT.rayTrace(1.0e9)
    thickAnalytic=pi*(RT.length*50)**2*bb(1.0e4,1.0e9)*1e23/(500*PC2CM)**2
#    assert almost_eq(thickim.sum()/1000,thickAnalytic, diff=0.1)
    print "optically thin test ok (ratio %.3f)"%(thickim.sum()/1000/thickAnalytic)

    thickPow=[log10(RT.rayTrace(x*1e9, suppressOutput=True).sum()) for x in xrange(1,11)]
    thickpf=polyfit([log10(x) for x in xrange(1,11)], thickPow,1)[0]
    assert almost_eq(thickpf, 2, 0.01)
    print "powerlaw spectrum for thick sphere",thickpf

    return RT

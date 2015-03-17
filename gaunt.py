from numpy import *
from utils import memoize, almost_eq
import math

def chebev(y,  c,  m):
    if (y < -1.0 or y > 1.0):     
        print ("x not in range")
        print ("returning 1")
        return 1
    d=0.
    dd=0.
    y2=2.*y
    for j in xrange(m-1,0,-1):
        sv=d
        d=y2*d-dd+c[j]
        dd=sv
    return y*d-dd+0.5*c[0]

def ghelp(t,z,nu):
    d = array([[8.986940175e00, -8.006936989e-1, -3.781305103e-1, 1.877213231e-2,  7.300158392e-2, -1.744671550e-3, -1.707268366e-2,  2.567331664e-4,  4.098322531e-3, 3.837562402e-5, -8.491991820e-4],        [-4.009515855e00,  9.466021705e-1,  1.102726322e-1, -1.004885705e-1,  3.576785497e-3,  2.864013856e-2, -4.694254776e-3, -9.155339970e-3,  1.635218463e-3, 2.938325230e-3, -3.615327726e-4],        [8.808871266e-1,  9.043402532e-2, -1.543619180e-2, -5.483366378e-2, -4.545307025e-3,  1.903394837e-2, 1.311691517e-3, -6.997479192e-3, -5.918883504e-4, 2.393747064e-3,  3.148015257e-4],        [2.640245111e-2, -9.608451450e-2,  8.310561114e-3, -4.520154409e-3, -1.017965604e-2,  7.091074494e-3, 5.316703136e-3, -3.571518641e-3, -2.333091048e-3, 1.328839809e-3,  8.909207650e-4],        [-4.580645915e-2, -1.885629865e-2,  2.179620525e-2, 8.366530426e-3, -9.530211924e-3, -9.668371391e-4,  5.178193095e-3, -2.096101038e-4, -2.484138313e-3,  9.135013312e-5,  9.869737522e-4],        [-3.568055702e-3,  1.050313890e-2,  4.259726289e-3,  3.700273930e-3, -3.450186162e-3, -2.999107465e-3,  2.451228935e-3,  1.553822487e-3, -1.359996060e-3, -7.137252303e-4,  6.134671184e-4],        [2.827798067e-3,  2.800889961e-3, -4.181588794e-3, 6.889320423e-4,  1.040482914e-3, -1.820642230e-3,-2.277321615e-5,  1.509584686e-3, -5.371426147e-5, 7.656848158e-4,  1.068883394e-4],        [3.365860195e-4, -1.078209202e-3, -1.770208330e-3, 9.460313195e-5,  1.407073544e-3, -3.874082085e-4, 8.182359057e-4,  6.212627837e-4,  5.553549563e-4,-3.504683798e-4, -2.046080100e-4]])
    gamma2 = 157833.0*z*z/t   # introduced gamma2, used a few times
    gamma = math.sqrt(gamma2)
    u = 4.797978e-11*nu/t
    HummerMask= (gamma2 >= 1.0e-3 and gamma2 <= 1.0e3 and u >= 1.0e-4 and u <= 31.622777)
    g=1
    if HummerMask:
        xu =math.log10(u)
        xu*=2.0
        xu+=2.5
        xu*=0.18181818181818181818
        
        xg = math.log10(gamma2)
        xg/=3.0
        
        cj =zeros(8,dtype=float)
        for j in xrange(8): cj[j]=(chebev(xg, d[j,:],11))
        g = chebev(xu,cj,8)
        

    # not in range for Hummer's fit - check to see if Scheuer's 
    # approximation is OK
    elif u < 1.0e-4 and gamma >= 1.0 :
    # use Scheuer's (1960) long-wavelength approximation (see Hummer)
    # this works ok for u < 10**-4 and gamma > 1
    #g[SchuMask]= -0.55133*(log(gamma[SchuMask]) + log(u[SchuMask]) + 0.056745)
        g=  math.log(gamma)
        g+= math.log(u)    #do individual ops to prevent uneccesary
        g+= 0.056745            #creation of large temporary arrays
        g*=-0.55133
    # not in range for Scheuer's fit, try Elwert's high-energy approx
    # (see Hummer) for u < 10**-4, and gamma < 1
    elif ((u < 1.0e-4) and  (gamma < 1.0)):
    # use Elwert's (1954) approximation (see Hummer)
        g=-math.log(u)
        g+=0.80888
        g*=0.55133
    #if none are applicable gaunt factor defaults to 1

    if g<0.1 : g=0.1 #floor gaunt factor at 0.1
        
    return g

gh=vectorize(memoize(lambda t,z,nu : ghelp(t,z,nu)))

def gaunt(t,z,nu):
    if almost_eq(t,t.flat[0], 1e-6).all():
        return ones_like(t)*ghelp(t.flat[0],z,nu)
    else:
        return gh(t,z,nu)

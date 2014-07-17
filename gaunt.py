from numpy import *
from scipy import constants as cns

def gaunt(t,z,nu):
    d = array([[8.986940175e00, -8.006936989e-1, -3.781305103e-1, 1.877213231e-2,  7.300158392e-2, -1.744671550e-3, -1.707268366e-2,  2.567331664e-4,  4.098322531e-3, 3.837562402e-5, -8.491991820e-4],        [-4.009515855e00,  9.466021705e-1,  1.102726322e-1, -1.004885705e-1,  3.576785497e-3,  2.864013856e-2, -4.694254776e-3, -9.155339970e-3,  1.635218463e-3, 2.938325230e-3, -3.615327726e-4],        [8.808871266e-1,  9.043402532e-2, -1.543619180e-2, -5.483366378e-2, -4.545307025e-3,  1.903394837e-2, 1.311691517e-3, -6.997479192e-3, -5.918883504e-4, 2.393747064e-3,  3.148015257e-4],        [2.640245111e-2, -9.608451450e-2,  8.310561114e-3, -4.520154409e-3, -1.017965604e-2,  7.091074494e-3, 5.316703136e-3, -3.571518641e-3, -2.333091048e-3, 1.328839809e-3,  8.909207650e-4],        [-4.580645915e-2, -1.885629865e-2,  2.179620525e-2, 8.366530426e-3, -9.530211924e-3, -9.668371391e-4,  5.178193095e-3, -2.096101038e-4, -2.484138313e-3,  9.135013312e-5,  9.869737522e-4],        [-3.568055702e-3,  1.050313890e-2,  4.259726289e-3,  3.700273930e-3, -3.450186162e-3, -2.999107465e-3,  2.451228935e-3,  1.553822487e-3, -1.359996060e-3, -7.137252303e-4,  6.134671184e-4],        [2.827798067e-3,  2.800889961e-3, -4.181588794e-3, 6.889320423e-4,  1.040482914e-3, -1.820642230e-3,-2.277321615e-5,  1.509584686e-3, -5.371426147e-5, 7.656848158e-4,  1.068883394e-4],        [3.365860195e-4, -1.078209202e-3, -1.770208330e-3, 9.460313195e-5,  1.407073544e-3, -3.874082085e-4, 8.182359057e-4,  6.212627837e-4,  5.553549563e-4,-3.504683798e-4, -2.046080100e-4]])
    gamma2 = 157833.0*z*z/t   # introduced gamma2, used a few times
    gamma = sqrt(gamma2)
    u = 4.797978e-11*nu/t
    HummerMask= logical_and(logical_and(logical_and(gamma2 >= 1.0e-3, gamma2 <= 1.0e3),u >= 1.0e-4),u <= 31.622777)
    g=ones_like(u)
    if HummerMask.any():
        xu = (2.0*log10(u)[HummerMask] + 2.5)*0.18181818181818181818
        xg = log10(gamma2[HummerMask])*0.333333333333333333333
        
        cj = [polyval(d[j,:],xg) for j in xrange(8)]
        g[HummerMask] = polyval(cj,xu)

    # not in range for Hummer's fit - check to see if Scheuer's 
    # approximation is OK
    SchuMask=logical_and(logical_not(HummerMask),logical_and(u < 1.0e-4, gamma >= 1.0))
    # use Scheuer's (1960) long-wavelength approximation (see Hummer)
    # this works ok for u < 10**-4 and gamma > 1
    g[SchuMask]= -0.55133*(log(gamma) + log(u) + 0.056745)[SchuMask]
    
    # not in range for Scheuer's fit, try Elwert's high-energy approx
    # (see Hummer) for u < 10**-4, and gamma < 1
    ElwertMask=logical_and(logical_not(logical_and(SchuMask,HummerMask)), logical_and((u < 1.0e-4), (gamma < 1.0)))
    # use Elwert's (1954) approximation (see Hummer)
    g[ElwertMask] = 0.55133*(0.80888 - log(u))[ElwertMask]
    #if none are applicable gaunt factor defaults to 1

    g[g<1]=1
        
    return g


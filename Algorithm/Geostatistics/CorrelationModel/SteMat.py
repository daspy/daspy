# h is the separation distance
# r is the range or distance parameter (r>0) which measures how quickly the correlations decay with distance
# v is the smoothness parameter (v>0)
# if v moves towards to infinity, the Matern function corresponds to a Gaussian model
# if v is equal to 0.5, the Matern function is Exponential
import numpy
from scipy.special import gamma,kv

def SteMat(h,r,kappa):

	# Matern isotropic covariance function (Stein's parameterization)

    n = numpy.size(h)
    corelation = numpy.zeros(n,dtype=numpy.double)
    
    for i in range(n):
        if h[i] == 0.0:
            corelation[i] = 1.0
        else:
            hr=numpy.double(h[i])/r
            
            bes=kv(kappa,2.0*numpy.sqrt(kappa)*hr)
            
            if not numpy.isfinite(bes):
                corelation[i] = 1.0
            elif bes == 0.0:
                corelation[i] = 0.0
            else:
                mult = 2.0**(1.0 - kappa)/gamma(kappa)*(2.0*numpy.sqrt(kappa)*hr)**kappa

                if not numpy.isfinite(mult):
                    corelation[i] = 0.0
                else:
                    corelation[i] = bes * mult
                    
    return corelation

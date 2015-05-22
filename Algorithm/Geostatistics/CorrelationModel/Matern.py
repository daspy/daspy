# h is the separation distance
# r is the range or distance parameter (r>0) which measures how quickly the correlations decay with distance
# kappa is the smoothness parameter (kappa>0)
# if kappa moves towards to infinity, the Matern function corresponds to a Gaussian model
# if kappa is equal to 0.5, the Matern function is Exponential
import numpy

from scipy.special import gamma,kv

def Matern(h,r,kappa):

	# Matern isotropic covariance function
    n = numpy.size(h)
    corelation = numpy.zeros(n,dtype=numpy.double)
    #if kappa == 0:
    #    corelation[:] = 0.0
    #else:
    for i in range(n):
        if h[i] == 0.0:
            corelation[i] = 1.0
        elif h[i] > 600 * r:
            corelation[i] = 0.0
        else:			
            hr = numpy.double(h[i])/r
            r1 = 2 ** (kappa-1.0) * gamma(kappa)
            bes = kv(kappa,hr)
            corelation[i] = (1.0/r1) * hr ** kappa * bes
    
    return corelation
    
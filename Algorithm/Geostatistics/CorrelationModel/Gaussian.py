# h is the separation distance
# r is the range or distance parameter (r>0) which measures how quickly the correlations decay with distance

import numpy

def Gaussian(h, r):

    n = numpy.size(h)
    corelation = numpy.zeros(n,dtype=numpy.double)
    for i in range(n):
        if h[i] == 0.0:
            corelation[i] = 1.0
        else:
            hr = numpy.double(h[i])/r
            corelation[i] = numpy.exp(-(hr * hr))

    return corelation

#print Gaussian([500000], [30000])
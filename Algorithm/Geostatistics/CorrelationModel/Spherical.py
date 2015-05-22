# h is the separation distance
# r is the range or distance parameter (r>0) which measures how quickly the correlations decay with distance

import numpy

def Spherical(h, r):
	
	n = numpy.size(h)
	corelation = numpy.zeros(n,dtype=numpy.double)
	for i in range(n):
		if h[i] == 0.0:
			corelation[i] = 1.0
		elif h[i] >= r:
			corelation[i] = 0.0
		else:
			hr = numpy.double(h[i])/r
			corelation[i] = 1.0 - hr * (1.5 - 0.5 * hr * hr)
			
	return corelation

#print Spherical([28000], 30000)

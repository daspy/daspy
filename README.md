The objective of DasPy development is for the mulisources and multivariate land data assimilation applications, such as soil moisture, soil temperature and joint state and parameter estimation. DasPy includes following components:

Model Operator:
CLM - Interface to Community Land Model 4.5
Assimilation Algorithm:
Local Ensemble Transform Kalman Filter (https://code.google.com/p/miyoshi/)
Observation Operator:
CMEM - Interface to Community Microwave Emission Modelling Platform 4.1 (https://software.ecmwf.int/wiki/display/LDAS/CMEM)
COSMIC - COsmic-ray Soil Moisture Interaction Code (http://cosmos.hwr.arizona.edu/)
TSF - Two-Source Formulation
Parallel Computing:
mpi4py - Message Passing Interface
parallelpython - Open Multi-Processing
Run Platform:
JUROPA (http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUROPA/JUROPA_node.html)
Linux
GCC (tested for 4.7.3, 4.8.2, 4.9.2)
Python 2.7
References:
Han, X., Franssen, H. J. H., Rosolem, R., Jin, R., Li, X., and Vereecken, H.: Correction of systematic model forcing bias of CLM using assimilation of cosmic-ray Neutrons and land surface temperature: a study in the Heihe Catchment, China, Hydrology and Earth System Sciences, 19, 615-629, 2015a.
Han, X., Li, X., Rigon, R., Jin, R., and Endrizzi, S.: Soil moisture estimation by assimilating L-band microwave brightness temperature with geostatistics and observation localization, Plos One, 10, e0116435, 2015b.
Han, X. J., Franssen, H. J. H., Montzka, C., and Vereecken, H.: Soil moisture and soil properties estimation in the Community Land Model with synthetic brightness temperature observations, Water Resour Res, 50, 6081-6105, 2014a.
Han, X. J., Jin, R., Li, X., and Wang, S. G.: Soil Moisture Estimation Using Cosmic-Ray Soil Moisture Sensing at Heterogeneous Farmland, Ieee Geoscience and Remote Sensing Letters, 11, 1659-1663, 2014b.
Han, X. J., Franssen, H. J. H., Li, X., Zhang, Y. L., Montzka, C., and Vereecken, H.: Joint Assimilation of Surface Temperature and L-Band Microwave Brightness Temperature in Land Data Assimilation, Vadose Zone J, 12, 0, 2013.
Han, X., Li, X., Franssen, H. J. H., Vereecken, H., and Montzka, C.: Spatial horizontal correlation characteristics in the land data assimilation of soil moisture, Hydrology and Earth System Sciences, 16, 1349-1363, 2012.
Montzka, C., Pauwels, V. R., Franssen, H. J., Han, X., and Vereecken, H.: Multivariate and multiscale data assimilation in terrestrial systems: a review, Sensors (Basel), 12, 16291-16333, 2012.
Han, X. and Li, X.: An evaluation of the nonlinear/non-Gaussian filters for the sequential data assimilation, Remote Sens Environ, 112, 1434-1449, 2008.
Acknowledgements:
The study of this work was supported by:

DFG (Deutsche Forschungsgemeinschaft) Forschergruppe 2131 "Data Assimilation for Improved Characterization of Fluxes across Compartmental Interfaces"
NSFC (National Science Foundation of China) project (grant number: 41271357, 91125001)
Transregional Collaborative Research Centre 32, financed by the German Science foundation
Supercomputing facilities of Forschungszentrum Julich (JUROPA)daspy

The objective of DasPy development is for the mulisources and multivariate land data assimilation applications, such as soil moisture, soil temperature and joint state and parameter estimation DasPy includes following components:

Model Operator:

•	CLM - Interface to Community Land Model 4.5

Assimilation Algorithm:

•	Local Ensemble Transform Kalman Filter (https://code.google.com/p/miyoshi/)

Observation Operator:

1.	CMEM - Interface to Community Microwave Emission Modelling Platform 4.1 (https://software.ecmwf.int/wiki/display/LDAS/CMEM)

2.	COSMIC - COsmic-ray Soil Moisture Interaction Code (http://cosmos.hwr.arizona.edu/)
3.	TSF - Two-Source Formulation
Parallel Computing:
1.	mpi4py - Message Passing Interface
2.	parallelpython - Open Multi-Processing
Run Platform:
1.	JUROPA (http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUROPA/JUROPA_node.html)
2.	Linux
3.	GCC (tested for 4.7.3, 4.8.2, 4.9.2)
4.	Python 2.7
References:
1.	Han, X., Franssen, H. J. H., Rosolem, R., Jin, R., Li, X., and Vereecken, H.: Correction of systematic model forcing bias of CLM using assimilation of cosmic-ray Neutrons and land surface temperature: a study in the Heihe Catchment, China, Hydrology and Earth System Sciences, 19, 615-629, 2015a.
2.	Han, X., Li, X., Rigon, R., Jin, R., and Endrizzi, S.: Soil moisture estimation by assimilating L-band microwave brightness temperature with geostatistics and observation localization, Plos One, 10, e0116435, 2015b.
3.	Han, X. J., Franssen, H. J. H., Montzka, C., and Vereecken, H.: Soil moisture and soil properties estimation in the Community Land Model with synthetic brightness temperature observations, Water Resour Res, 50, 6081-6105, 2014a.
4.	Han, X. J., Jin, R., Li, X., and Wang, S. G.: Soil Moisture Estimation Using Cosmic-Ray Soil Moisture Sensing at Heterogeneous Farmland, Ieee Geoscience and Remote Sensing Letters, 11, 1659-1663, 2014b.
5.	Han, X. J., Franssen, H. J. H., Li, X., Zhang, Y. L., Montzka, C., and Vereecken, H.: Joint Assimilation of Surface Temperature and L-Band Microwave Brightness Temperature in Land Data Assimilation, Vadose Zone J, 12, 0, 2013.
6.	Han, X., Li, X., Franssen, H. J. H., Vereecken, H., and Montzka, C.: Spatial horizontal correlation characteristics in the land data assimilation of soil moisture, Hydrology and Earth System Sciences, 16, 1349-1363, 2012.
7.	Montzka, C., Pauwels, V. R., Franssen, H. J., Han, X., and Vereecken, H.: Multivariate and multiscale data assimilation in terrestrial systems: a review, Sensors (Basel), 12, 16291-16333, 2012.
8.	Han, X. and Li, X.: An evaluation of the nonlinear/non-Gaussian filters for the sequential data assimilation, Remote Sens Environ, 112, 1434-1449, 2008.
Acknowledgements:
The study of this work was supported by:
1.	DFG (Deutsche Forschungsgemeinschaft) Forschergruppe 2131 "Data Assimilation for Improved Characterization of Fluxes across Compartmental Interfaces"
2.	NSFC (National Science Foundation of China) project (grant number: 41271357, 91125001)
3.	Transregional Collaborative Research Centre 32, financed by the German Science foundation
4.	Supercomputing facilities of Forschungszentrum Julich (JUROPA)


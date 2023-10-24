This folder contains the following Matlab scripts and data files:

#

- 'main_realizations_implicitFDM.m' and 'main_realizations_explicitFDM.m'
are codes to compute velocity realizations with implicit/explicit FDM, respectively.

- 'implicitFDM.m' is the version for homogeneous boundary value problems of the FDM solver published in 
https://github.com/CDAlecsa/Boundary-Value-Problem-2D.

- 'explicitFDM.m' is the "deterministic GRW" solver 'realiz_Gauss_GRW.m' for advection-diffusion transport from 
https://github.com/PMFlow/RichardsEquation/tree/main/2D/Saturated_2D_FlowTransport.

- 'Kraichnan_Gauss_param.m' and 'K' are versions of the functions used to compute random phases and wavenumbers and realizations of the hydraulic conductivity, respectively, from
https://github.com/PMFlow/FlowBenchmark/tree/master/K-field_generators/Matlab.

- 'initstate.m' fixes the seed of the random number generator 'randn' used in 'Kraichnan_Gauss_param.m'.

- 'test_implicitFDM' and 'test_explicitFDM' are folders containing ensembles of 100 realizations of a random velocity field computed with the implicit and explicit FDM, respectively.



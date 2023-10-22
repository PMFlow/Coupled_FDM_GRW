
Solutions of stationary flows in saturated porous media with hydraulic conductivity modeled as a spatial random function are computed with numerical schemes based on implicit and explicit finite difference methods (FDM) designed in (https://doi.org/10.1016/j.advwatres.2020.103558, Section 3.1 and 3.6).

Advection-diffusion solutions of passive transport in aquifers are otained with the unbiased global random walk algorithm (GRW) (https://doi.org/10.1016/j.advwatres.2021.103935, Section 4.2.2).

Monte Carlo inferences of the velocity of the center of mass of the solute plume and of the ensemble dispersion coefficients are obtained from GRW solutions computed with realizations of velocity fields from ensembles of FDM flow solutions.

The coupled FDM-GRW approach using the implicit FDM scheme is found to be at least two orders of magnitude faster than the approach based on the explicit FDM scheme. The results obtained by both approaches are close to each other and close to the results of the linear approximation using the Kraichnan generator of the velocity field. Note that the linear approximation is one order of magnitude faster than the implicit FDM but it is limited to small variances of the lnK field.

The Matlab codes and data are contained in the following folders:

FDM_Flow
- implicit and explicit FDM flow solutions and ensembles of realizations.

GRW_Transport
- the computation of the center of mass velcoity and of the dispersion coefficients.

Comparisons
-comparisons between FDM approaches and a linear approximation.


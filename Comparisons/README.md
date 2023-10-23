This folder contains the following Matlab scripts and data files:

#

- 'velocity_samples.m' compares realizations of the velocity field computed with the implicit and explicit FDM, respectively.

- 'comparison_coeff_expl_impl.m' is the Matlab script that compares center of mass velocities and dispersion coefficients inferred from the ensembles of 100 realizations of the velocity field from

 '..\GRW_Transport\ensemble_coefficients_implicitFEM' and 

 '..\GRW_Transport\ensemble_coefficients_explicitFEM'.

- 'comparison_coeff_FDM_linear.m' compares center of mass velocities and dispersion coefficients inferred from ensembles of 10^4 realizations of the advection-diffusion transport simulated with velocity fields computed with implicit FDM to the same quantities inferred from 10^4 realizations of the Kraichnan generator.

- 'ensemble_coefficients_explicitFD_Nmod100R1e4.mat' is the data file containing coefficients based on implicit FDM.

- 'linear_approximation.mat' is the data file containing coefficients derived with the Kraichnan procedure.

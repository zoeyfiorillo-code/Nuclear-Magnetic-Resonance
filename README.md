# Nuclear Magnetic Resonance Data Analysis
Two programs were written in MATLAB to determine the spin-lattice ($T_1$), spin-spin ($T_2$), and effective transverse ($T_2^*$) relaxation times from experimental nuclear magnetic resonance (NMR) data. 

## NMR.m
The objective of the program was to determine the spin-lattice ($T_1$), spin-spin ($T_2$), and effective transverse ($T_2^*$) relaxation times. It outputs three graphs containing the original data, fitted curves, equations with fitted parameters, and coefficients of determination ($R^2$). It also prints the fitted parameters and their uncertainties. 

## NMR_Fit_Comp.m
The objective of the program was to compare how nonlinear exponential curve fitting and linearized regression determined the effective transverse ($T_2^*$) relaxation time. It outputs two graphs containing the original data, fitted curves, equations with fitted parameters, and chi-squared ($\chi^2$) statistics. It also prints the fitted parameters and their uncertainties. 

## Dependencies and Acknowledgements
The linear regression analysis requires an additional MATLAB program, linreg.m, developed by Alejandro L. Garcia (Copyright © 1999) for *Numerical Methods for Physics*, 2nd Edition. 

The file is available separately through the original GitHub repository: [NM4P](https://github.com/AlejGarcia/NM4P/blob/5d071c8044f633602c2eff60487f91918c97170c/MatlabRevised/linreg.m)

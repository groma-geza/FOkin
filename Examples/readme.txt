The following examples illustrate all the major capabilities of the FOkin toolbox, by applying the algorithms and reproducing the figures and tables presented in [1]. Note, that due to the stochastic nature of cross-validation and Bayesian optimization, the values of the calculated data can be slightly different for every run.

example1.m – script applying Algorithm2 (based on the RCV(nv) version of cross-validation) for the analysis of simulated data derived from a complex model of the bacteriorhodopsin photocycle.

example2.m – script for excluding distributed kinetics on the data analyzed by example1.m by model selection based on 10-fold CV (first step of Algorithm 3).

example3.m – script applying Algorithm2 for the analysis of experimental ultrafast fluorescence kinetic data measured on the coenzyme FAD.

example4.m - script for excluding distributed kinetics on the data analyzed by example3.m.

example5.m – script for analysis of simulated data with distributed kinetics by both RCV(nv) and 10-fold CV.

example6.m – script demonstrating the differences in the results of model selections executed without cross-validation, with 10-fold CV and with RCV(nv).

example7.m – script for analysis of simulated data with realistic noise.

example8.m – script for analysis of simulated data of Erlang distribution without and with exponential components.

example9.m – script for analysis of simulated data of second-order kinetics without and with exponential components.


1. Zimányi L, Sipos Á, Sarlós F, Nagypál R, Groma GI. Machine-learning model selection and parameter estimation from kinetic data of complex first-order reaction systems. PLoS One. 2021.
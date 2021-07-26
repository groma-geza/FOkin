clear
fname = 'example4.mat';

%% Load result of example3.m
load('example3.mat') % all data (including t0 and fwhm) are stored in a 
                     % FOkin object
%% Execute 10-fold CV in the 2D space of lambda end omega

lambda = [1.E-6, 1.E-1];
omega = [1.E-8, 1];
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', fokin.name, '''...'])
selected_dual = fokin.select_model(lambda, omega, maxiter, info, 'kCV');
selected_dual.do_sections(.4, .5, 1)
selected_dual.show()

save(fname, 'fokin', 'selected_dual')

%% Recalculate 10-fold CV in the omega space with the optimal value of lambda.
load(fname)
lambda = selected_dual.best_lambda; % optimal value of lambda
omega = [1.E-8, 1];
maxiter = 400;
info = selected_dual.info;
selected_omega = fokin.select_model(lambda, omega, maxiter, info, 'kCV');
selected_omega.show()
save(fname, 'fokin', 'selected_omega', '-append')
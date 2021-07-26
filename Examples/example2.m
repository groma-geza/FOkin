clear
fname = 'example2.mat';

%% Load simulated data and corrupt it by noise
noise_sigma = 1.E-3; % noise level
wavelength_step = 10; % wavelength step

[data, time, wavelength] = create_bR_data(noise_sigma, wavelength_step, 0);

%% Set up a FOkin object
name = ['  sigma = ', num2str(noise_sigma, '%.1E')];            
fokin = FOkin(data, time, wavelength);
fokin.name = name;

fopt =  fokin.options;
fopt.num_par_workers = Inf;
fopt.extension_hi = 0;
fopt.extension_lo = 0.1;
fopt.signal_label = '\DeltaA (rel)';

optimizer = fokin.optimizer;
eopt = optimizer.options;
eopt.MAX_ITER = 1.E6;
eopt.ABSTOL   = 1E-11;
eopt.RELTOL   = 1E-7;

%% Execute 10-fold CV in the 2D space of lambda end omega
% lambda = [1.E-7 1.E-2]; % for sigma = 1.E-4
lambda = [5.E-6, 1.E-1]; % for sigma = 1.E-3
% lambda = [5.E-5, 2.E-1]; % for sigma = 1.E-2
omega = [1.E-7, 1];
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', name, '''...'])
selected_dual = fokin.select_model(lambda, omega, maxiter, info, 'kCV');
selected_dual.show()
selected_dual.do_sections(.4, .4, 1) % thin slices for visualization in 1D
save(fname, 'fokin', 'selected_dual')

%% Recalculate 10-fold CV in the omega space with the optimal value of lambda
load(fname)
lambda = selected_dual.best_lambda; % optimal value of lambda
omega = [1.E-7, 1];
maxiter = 100;
info = selected_dual.info;
selected_omega = fokin.select_model(lambda, omega, maxiter, info, 'kCV');
selected_omega.show()
save(fname, 'fokin', 'selected_omega', '-append')
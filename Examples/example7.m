clear
fname = 'example7.mat';

%% Load simulated data and corrupt it by noise
wavelength_step = 10; % wavelength step

[data, time, wavelength, aver_noise_level] = ...
    create_bR_data_with_real_noise(wavelength_step, 1);

%% Set up a FOkin object
name = '  real noise';            
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

%% Execute RCVnv in the 2D space of lambda end omega
lambda = [1.E-5, 1];
omega = [1.E-7, 1];
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', name, '''...'])
selected_dual = fokin.select_model(lambda, omega, maxiter, info, 'RCVnv');
selected_dual.show()

selected_dual.do_sections(.4, .4, 1) % thin slices for visualization in 1D
save(fname, 'fokin', 'selected_dual')

%% Recalculate RCVnv in the omega space with the optimal value of lambda
load(fname)
lambda = selected_dual.best_lambda; % optimal value of lambda
omega = [1.E-7, 1];
maxiter = 100;
info = selected_dual.info;
selected_omega = fokin.select_model(lambda, omega, maxiter, info, 'RCVnv');
selected_omega.count_support(100, 1, 1) % support size for visualization
selected_omega.show()
save(fname, 'fokin', 'selected_omega', '-append')

%% Recalculate RCVnv in the lambda space with the optimal value of omega
load(fname)
lambda = [1.E-05 1];
omega = selected_omega.best_omega; % optimal value of omega
maxiter = 100;
info = selected_omega.info;
selected_lambda = fokin.select_model(lambda, omega, maxiter, info,...
    'RCVnv');
selected_lambda.count_features(100, 1) % number of features
                                       % for visualization
selected_lambda.show()
save(fname, 'fokin', 'selected_lambda', '-append')

%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda.best_lambda;
omega = selected_omega.best_omega;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = omega;
info = selected_lambda.info;
fokin.optimize()
result = fokin.results(info);
result.denorm()
result.show(416)
discrete = result.discretize();
discrete.showd() % S3 Table
save(fname, 'fokin', 'discrete', '-append')



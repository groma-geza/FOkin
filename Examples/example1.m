clear
fname = 'example1.mat';

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

%% Execute RCVnv in the 2D space of lambda end omega
% lambda = [5.E-7, 1.E-1]; % for sigma = 1.E-7
% lambda = [1.E-6, 1.E-1]; % for sigma = 1.E-6
% lambda = [1.E-6, 1.E-1]; % for sigma = 1.E-5
% lambda = [1.E-5, 1.E-1]; % for sigma = 1.E-4
lambda = [1.E-5, 1]; % for sigma = 1.E-3
% lambda = [1.E-4, 1.E1]; % for sigma = 1.E-2

omega = [1.E-7, 1];
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', name, '''...'])
selected_dual = fokin.select_model(lambda, omega, maxiter, info, 'RCVnv');
selected_dual.show() % Fig 3A

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
selected_omega.show() % Fig 3B
save(fname, 'fokin', 'selected_omega', '-append')

%% Recalculate RCVnv in the lambda space with the optimal value of omega
load(fname)
% lambda = [4.E-07 2.E-5]; % for sigma = 1.E-7
% lambda = [7.E-06 1.5E-4]; % for sigma = 1.E-6
% lambda = [1.E-07 1.E-1]; % for sigma = 1.E-5
% lambda = [2.E-04 5.E-2]; % for sigma = 1.E-4 % very narrow valley to find
lambda = [1.E-06 1]; % for sigma = 1.E-3
% lambda = [1.E-04 1.E1]; % for sigma = 1.E-2
omega = selected_omega.best_omega; % optimal value of omega
maxiter = 100;
info = selected_omega.info;
selected_lambda = fokin.select_model(lambda, omega, maxiter, info,...
    'RCVnv');
selected_lambda.count_features(100, 1) % number of features
                                       % for visualization
selected_lambda.show() % Fig 3C
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
result.show(416) % Fig 4A, S4 Fig A-C
discrete = result.discretize();
discrete.showd()
discrete.neglect(0.05)
discrete.showd() % Fig 4B, Table 1, S2 Table
save(fname, 'fokin', 'discrete', '-append')
 
%% Execute exponential fitting
load(fname)
discrete.do_expfit()
discrete.showexp(416) % Fig 4C, Table 1, S2 Table, S4 Fig D, 
                      % Fig 4A blue arrows
save(fname, 'fokin', 'discrete', '-append')

discrete.reset()
discrete.do_expfit()
discrete.showexp(416) % S4 Fig E, S5 Fig
 


clear
fname = 'example5.mat';

%% Load simulated data and corrupt it by noise
sigma = 1.E-3;

[time, data] = create_distributed_data(sigma, 0);

%% Set up a FOkin object
group_param = 0; %dummy value (the simulated data fall into a single group) 
fokin = FOkin(data, time, group_param);
fokin.name = ['distributed  sigma = ', num2str(sigma, '%.1E')];

fopt = fokin.options;
fopt.extension_hi = 0;
fopt.extension_lo = 0;
fopt.num_par_workers = Inf;
fopt.group_lasso = 0;
fopt.signal_label = 'Signal';
fopt.time_unit = 'a.u.';

% In create_distributed_data k = exp(-E/50). 
% The corresponding n_tau_decade is:
fopt.n_tau_decade = 50/log10(exp(1)); 
fopt.nonzero_limit = 1.E-6;
fopt.min_thread_length = 1;

optimizer = MultiElnetPDCO;
eopt = optimizer.options;
eopt.OptTol = 1.E-15;
eopt.MaxIter = 1.E5;

fokin.optimizer = optimizer;

%% Execute RCVnv in the 2D space of lambda end omega
% lambda = [1.E-3, 2.E0]; % for sigma = 1.E-2
lambda = [1.E-4, 1.E1]; % for sigma = 1.E-3
% lambda = [1.E-4, 1.E1]; % for sigma = 1.E-4
omega = [1.E-7, 1 - 1.E-6]; % PDCO does not converge well with omega = 1
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', fokin.name, '''...'])
selected_dual_RCVnv = fokin.select_model(lambda, omega, maxiter, info,...
    'RCVnv');
selected_dual_RCVnv.show() % Fig 9A left
selected_dual_RCVnv.do_sections(.3, .3, 1) % thin slices for visualization 
                                           %  in 1D
save(fname, 'fokin', 'selected_dual_RCVnv')

%% Recalculate RCVnv in the omega space with the optimal value of lambda
load(fname)
lambda = selected_dual_RCVnv.best_lambda; % optimal value of lambda
omega = [1.E-7, 1 - 1.E-6];
maxiter = 100;
info = selected_dual_RCVnv.info;
selected_omega_RCVnv = fokin.select_model(lambda, omega, maxiter, info,...
    'RCVnv');
selected_omega_RCVnv.count_support(100, 0, 1) % support size for 
                                        % visualization
                                        % Does not include zero since PDCO
                                        % requires nonzero L2 penalty.
selected_omega_RCVnv.show() % Fig 9B left
save(fname, 'fokin', 'selected_omega_RCVnv', '-append')

%% Recalculate RCVnv in the lambda space with the optimal value of omega
load(fname)
% lambda = [1.E-3, 2.E0]; % for sigma = 1.E-2
lambda = [1.E-4, 1.E1]; % for sigma = 1.E-3
% lambda = [1.E-4, 1.E1]; % for sigma = 1.E-4
omega = selected_omega_RCVnv.best_omega; % optimal value of omega
maxiter = 100;
info = selected_omega_RCVnv.info;
selected_lambda_RCVnv = fokin.select_model(lambda, omega, maxiter, info,...
    'RCVnv');
selected_lambda_RCVnv.count_features(100, 1) % number of features
                                             % for visualization
selected_lambda_RCVnv.show() % Fig 9C left 
save(fname, 'fokin', 'selected_lambda_RCVnv', '-append' )

%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda_RCVnv.best_lambda;
omega = selected_omega_RCVnv.best_omega;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = omega;
info = selected_lambda_RCVnv.info;
fokin.optimize()
result_RCVnv = fokin.results(info);
result_RCVnv.show() % Fig 10A-B left
save(fname, 'fokin', 'result_RCVnv', '-append')

%% Execute 10-fold CV in the 2D space of lambda end omega.
load(fname)
% lambda = [1.E-3, 2.E0]; % for sigma = 1.E-2
lambda = [1.E-5, 1.E-1]; % for sigma = 1.E-3
% lambda = [1.E-6, 1.E-2]; % for sigma = 1.E-4
omega = [1.E-7, 1];
maxiter = 600;
info.dummy = 'dummy';

disp(['Running model selection for ''', fokin.name, '''...'])
selected_dual_kCV = fokin.select_model(lambda, omega, maxiter, info,...
    'kCV');
selected_dual_kCV.show() % Fig 9A right
selected_dual_kCV.do_sections(.3, .3, 1)
save(fname, 'fokin',  'selected_dual_kCV', '-append')

%% Recalculate 10-fold CV in the omega space with the optimal value of lambda.
load(fname)
lambda = selected_dual_kCV.best_lambda;
omega = [1.E-7, 1];
maxiter = 300;
info = selected_dual_kCV.info;
selected_omega_kCV = fokin.select_model(lambda, omega, maxiter, info,...
    'kCV');
selected_omega_kCV.count_support(100, 0, 1)
selected_omega_kCV.show() % Fig 9B right
save(fname, 'fokin', 'selected_omega_kCV', '-append')

%% Recalculate 10-fold CV in the lambda space with the optimal value of omega.
load(fname)
% lambda = [1.E-3, 2.E0]; % for sigma = 1.E-2
lambda = [1.E-5, 1.E-1]; % for sigma = 1.E-3
% lambda = [1.E-6, 1.E-2]; % for sigma = 1.E-4
omega = selected_omega_kCV.best_omega;
maxiter = 100;
info = selected_omega_kCV.info;
selected_lambda_kCV = fokin.select_model(lambda, omega, maxiter, info,...
    'kCV');
selected_lambda_kCV.count_features(100, 1)
selected_lambda_kCV.show() % Fig 9C right
save(fname, 'fokin', 'selected_lambda_kCV', '-append')

%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda_kCV.best_lambda;
omega = selected_omega_kCV.best_omega;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = omega;
info = selected_lambda_kCV.info;
fokin.optimize()
result_kCV = fokin.results(info);
result_kCV.show() % Fig 10A-B right, Fig 11 blue, gold an brow lines
save(fname, 'fokin', 'result_kCV', '-append')

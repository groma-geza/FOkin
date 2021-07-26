clear
fname = 'example9.mat';

C = 8;
noise_sigma = 1.E-7;
with_exp = 1;
show = 0;

[time, data, txt1] = create_2nd_order_data(C, noise_sigma, with_exp, show);

%% Set up a FOkin object
group_param = 0; %dummy value (the simulated data fall into a single group) 
fokin = FOkin(data, time, group_param);
name = [txt1, '  sigma = ', num2str(noise_sigma, '%.1E')];
fokin.name = name;

fopt = fokin.options;
fopt.extension_hi = 0;
fopt.extension_lo = 0;
fopt.num_par_workers = Inf;
fopt.group_lasso = 0;
fopt.signal_label = 'Signal';
fopt.time_unit = 'a.u.';

fopt.min_thread_length = 1;

optimizer = MultiElnetADMM;
eopt = optimizer.options;
eopt.MAX_ITER = 1.E6;
eopt.ABSTOL   = 1E-11;
eopt.RELTOL   = 1E-7;

fokin.optimizer = optimizer;

if isfile(fname)
    save(fname, 'fokin', '-append')
else
    save(fname, 'fokin')
end

%% Execute RCVnv in the 2D space of lambda end omega
load(fname)
lambda = [1.E-8, 1]; % for sigma = 1.E-3 with two exponential components
omega = [1.E-8, 1];
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', fokin.name, '''...'])
selected_dual = fokin.select_model(lambda, omega, maxiter, info, 'RCVnv');
selected_dual.show()

selected_dual.do_sections(.6, .6, 1) % thin slices for visualization in 1D
save(fname, 'fokin', 'selected_dual', '-append')

%% Recalculate RCVnv in the omega space with the optimal value of lambda
load(fname)
lambda = selected_dual.best_lambda; % optimal value of lambda
omega = [1.E-8, 1];
maxiter = 100;
info = selected_dual.info;
selected_omega = fokin.select_model(lambda, omega, maxiter, info, 'RCVnv');
selected_omega.count_support(100, 1, 1) % support size for visualization
selected_omega.show()
save(fname, 'fokin', 'selected_omega', '-append')

%% Recalculate RCVnv in the lambda space with the optimal value of omega
load(fname)
lambda = [1.E-8, 1]; 
omega = selected_omega.best_omega; % optimal value of omega
maxiter = 100;
info = selected_omega.info;
selected_lambda = fokin.select_model(lambda, omega, maxiter, info,...
    'RCVnv');
selected_lambda.count_features(100, 1) % number of features
                                       % for visualization
selected_lambda.show()
save(fname, 'fokin', 'selected_lambda', '-append')
% 
%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda.best_lambda;
omega = selected_omega.best_omega;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = omega;
info.dummy = 'dummy' ;
fokin.optimize()
result = fokin.results(info);
result.show() % Fig 13 A or C, S8 Fig A or C
save(fname, 'fokin', 'result', '-append')

%% Execute 10-fold CV in the 2D space of lambda end omega
load(fname)
lambda = [1.E-12, 1]; % sigma 1.E-7 witht exp
% lambda = [1.E-8, 1]; % everything else
omega = [1.E-8, 1];
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', fokin.name, '''...'])
selected_dual_kcv = fokin.select_model(lambda, omega, maxiter, info, 'kCV');
selected_dual_kcv.show()

selected_dual_kcv.do_sections(.6, .6, 1) % thin slices for visualization in 1D
save(fname, 'fokin', 'selected_dual_kcv', '-append')

%% Recalculate 10-fold CV in the omega space with the optimal value of lambda
load(fname)
lambda = selected_dual_kcv.best_lambda; % optimal value of lambda
omega = [1.E-8, 1];
maxiter = 300;
info = selected_dual_kcv.info;
selected_omega_kcv = fokin.select_model(lambda, omega, maxiter, info, 'kCV');
% selected_omega_kcv.count_support(100, 1, 1) % support size for visualization
selected_omega_kcv.show()
save(fname, 'fokin', 'selected_omega_kcv', '-append')

%% Recalculate 10-fold CV in the lambda space with the optimal value of omega
load(fname)
lambda = [1.E-12, 1]; % sigma 1.E-7 witht exp
% lambda = [1.E-8, 1]; % everything else
omega = selected_omega_kcv.best_omega; % optimal value of omega
maxiter = 100;
info = selected_omega_kcv.info;
selected_lambda_kcv = fokin.select_model(lambda, omega, maxiter, info,...
    'kCV');
% selected_lambda_kcv.count_features(100, 1) % number of features
                                       % for visualization
selected_lambda_kcv.show()
save(fname, 'fokin', 'selected_lambda_kcv', '-append')
% 
%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda_kcv.best_lambda;
omega = selected_omega_kcv.best_omega;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = omega;
info.dummy = 'dummy' ;
fokin.optimize()
result_kcv = fokin.results(info);
result_kcv.show() % Fig 12 B or D, S8 Fig B or D
save(fname, 'fokin', 'result_kcv', '-append')

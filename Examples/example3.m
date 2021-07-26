clear
fname = 'example3.mat';

%% Load data
load ('FAD_data.mat');
name = 'FAD fluorescence';
data_norm = data/max(max(data));
figure, semilogx(delay, data_norm) % Fig 5A
xlabel('Delay (ps)')
y_lab = 'Fluorescence (rel)';
ylabel(y_lab)
title(name)
figure, plot(wavelength, data_norm) % Fig 5B
wav_lab = 'Wavelength (nm)';
xlabel(wav_lab)
ylabel(y_lab)
title(name)

%% Determine instrument response function parameters
ropt = FOkin.RisefitOptions;
ropt.separate_norm = true;
ropt.elnet_lambda = 5.E-2;
ropt.maxiter = 300;
ropt.show_details = 3;
ropt.title_text = name;

fopt = FOkinOptions;
fopt.extension_lo = 0;
fopt.num_par_workers = Inf;
fopt.n_t0_knots = 4;
fopt.t0_range = [0.5, 1.2];
fopt.fwhm_range = [0.2, 0.4];

disp(['Running risefit for ''', name, '''...'])
[t0, fwhm] = FOkin.risefit(data, delay, wavelength, ropt, fopt);

%% Create a FOkin object
fokin = FOkin(data, delay, wavelength);
fokin.t0 = t0;
fokin.fwhm = fwhm;
fokin.name = name;
fokin.options = fopt;
save(fname, 'fokin')

%% Execute RCVnv in the 2D space of lambda end omega
load(fname)
lambda = [1.E-4, 1.E1];
omega = [1.E-8, 1];
maxiter = 400;
info.dummy = 'dummy'; % dummy input info

rng('shuffle') % shuffle random number generator

disp(['Running model selection for ''', fokin.name, '''...'])
selected_dual = fokin.select_model(lambda, omega, maxiter, info);
selected_dual.show() % Fig 6A

selected_dual.do_sections(.4, .5, 1) % thin slices for visualization in 1D
save(fname, 'fokin', 'selected_dual', '-append')

%% Recalculate RCVnv in the omega space with the optimal value of lambda
load(fname)
lambda = selected_dual.best_lambda; % optimal value of lambda
omega = [1.E-8, 1];
maxiter = 100;
info = selected_dual.info;
selected_omega = fokin.select_model(lambda, omega, maxiter, info);
selected_omega.count_support(100, 1, 1) % support size for visualization
selected_omega.show() % Fig 6B
save(fname, 'fokin', 'selected_omega', '-append')

%% Recalculate RCVnv in the lambda space with the optimal value of omega
load(fname)
lambda = [2.E-4, 2];
omega = selected_omega.best_omega; % optimal value of omega
maxiter = 100;
info = selected_omega.info;
selected_lambda = fokin.select_model(lambda, omega, maxiter, info);
selected_lambda.count_features(100, 1) % number of features
                                       % for visualization
selected_lambda.show() % Fig 6C
save(fname, 'fokin', 'selected_lambda', '-append')

%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda.best_lambda;
omega = selected_omega.best_omega;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = omega;
info = selected_lambda.info;
disp(['Estimating kinetic parameters for ''', fokin.name, '''...'])
fokin.optimize()
result = fokin.results(info);
result.show(520) % Fig 7A, S6 Fig A-C

discrete = result.discretize();
discrete.neglect(5.E-3)
discrete.options.smooth_val = 1.E-2; % deep copy of fokin.options
discrete.showd() % Fig 7B, Table 2
save(fname, 'fokin', 'discrete', '-append')

%% Execute exponential fitting
load(fname)
disp(['Executing exponential fitting for ''', fokin.name, '''...'])
discrete.do_expfit()
discrete.options.thread_colors =...
    [0.39 0.83 0.07   % green    
    0   0    1      % blue
    1    0    0      % red
    0.06 1    1      % cyan
    1    0    1      % magenta
    0.93 0.69 0.13   % gold
    0.64 0.08 0.18]; % brown
discrete.showexp(520) % Fig 7C, Table 2, Fig 7A blue arrows, S6 Fig D
save(fname, 'fokin', 'discrete', '-append')
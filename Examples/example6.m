clear
fname = 'example6.mat';

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

%% Calculate MSE on a grid of lambda with omega = 0 without CV
optimizer.omega = 0;
lambda_range = [-6, -1];
lambda_grid = 10.^(lambda_range(1):.1:lambda_range(2));

info.sigma = noise_sigma;
info.lambda_grid = lambda_grid;
info.omega = 0;

n = length(lambda_grid);
MSE_total = NaN(1, n);
hfig = figure; % Fig 2A cyan curve

for ii = 1:n
    optimizer.lambda = lambda_grid(ii);
    fokin.optimize();
    res = fokin.results(info);
    MSE_total(ii) = res.MSE_total;
    loglog(lambda_grid, MSE_total, 'c.-', 'LineWidth', 1, 'MarkerSize', 10)
    xlabel('\lambda')
    ylabel('MSE')
    drawnow
end
save(fname, 'fokin', 'lambda_grid', 'MSE_total', 'info')

%% Calculate 10-fold CV in the lambda space with omega =0
load(fname)
lambda = [1.E-06 1.E-1];
omega = 0;
maxiter = 100;

rng('shuffle') % shuffle random number generator

selected_lambda1 = fokin.select_model(lambda, omega, maxiter, info, 'kCV');
selected_lambda1.show() % Fig 2A blue and red curves
save(fname, 'fokin', 'selected_lambda1', '-append')

%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda1.best_lambda;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = 0;
info = selected_lambda1.info;
fokin.optimize()
result1 = fokin.results(info);
result1.show(416) % Fig 2B
save(fname, 'fokin', 'result1', '-append')

%% Calculate RCVnv in the lambda space with omega =0
load(fname)
lambda = [1.E-06 1.E0];
omega = 0;
maxiter = 100;
selected_lambda2 = fokin.select_model(lambda, omega, maxiter, info,...
    'RCVnv');
selected_lambda2.show() % Fig 2C
save(fname, 'fokin', 'selected_lambda2', '-append')

%% Estimate the kinetic parameters
load(fname)
lambda = selected_lambda2.best_lambda;
optimizer = fokin.optimizer;
optimizer.lambda = lambda;
optimizer.omega = 0;
info = selected_lambda2.info;
fokin.optimize()
result2 = fokin.results(info);
result2.show(416) % Fig 2D
save(fname, 'fokin', 'result2', '-append')
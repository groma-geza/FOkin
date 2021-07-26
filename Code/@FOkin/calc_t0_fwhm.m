function [t0, fwhm, rel_params] = calc_t0_fwhm(self, maxiter)
arguments
    self
    maxiter(1,1) double {mustBeInteger, mustBePositive}
end

data = self.data;
time = self.time;
x = self.group_param;
p = self.p;

% Create a deep copy of the options to make the following changes temporal.
opt = self.options.copy(); 
opt.bo_IsObjectiveDeterministic = true;
if opt.num_par_workers
    opt.bo_UseParallel = true;
else
    opt.bo_UseParallel = false;
end

% Get Name, Value input parameter pairs for bayesopt() from the properties
% of opt.
controls = get_BO_controls(opt, maxiter);

n_t0_knots = opt.n_t0_knots;
err = false;
if (floor(n_t0_knots) == n_t0_knots) % integer
    if ((n_t0_knots < 1) || (n_t0_knots > p))
        err = true;
    end
else
    err = true;
end
if err
    error(['Possile values of ''n_t0_knots'' are positive integers ',...
        '<= ''self.p''.'])
end

n_fwhm_knots = opt.n_fwhm_knots;
err = false;
if (floor(n_fwhm_knots) == n_fwhm_knots) % iteger
    if ((n_fwhm_knots < 1) || (n_fwhm_knots > p))
        err = true;
    end
else
    err = true;
end
if err
    error(['Possile values of ''n_fwhm_knots'' are positive integers ',...
        '<= ''self.p''.'])
end

t0_range = opt.t0_range;
fwhm_range = opt.fwhm_range;
fwhm_fixed = opt.fwhm_fixed;

num_par_workers = opt.num_par_workers;
if ((num_par_workers > 0) && isempty(gcp('nocreate')))
    if isinf(num_par_workers)
        parpool(opt.parpool_spec);
    else
       parpool(opt.parpool_spec, num_par_workers); 
    end
end

% Estimate t0.
t0_raw = zeros(1, p);
for ii = 1:p
    data_ii = data(:,ii);
    [maxval, idx] =  max(abs(data_ii));
    data_ii = data_ii(1:idx)/maxval;
    time_ii = time(1:idx);
    lo = log10(time_ii(1));
    hi = log10(time_ii(end));
    step = (hi - lo)/999;
    fine_time_ii = 10.^(lo:step:hi);
    fine_data_ii = interp1(time_ii, data_ii, fine_time_ii, 'pchip');
    level = 0.75; % Mean of 1 and 0.5 corresonding to decay time constant
                  % of 0 and Inf, respectively
    idx = find(abs(fine_data_ii) < level, 1, 'last');
    t0_raw(ii) = mean([fine_time_ii(idx), fine_time_ii(idx + 1)]);
end

t0_knot_idx = get_knot_indices(p, n_t0_knots);
t0_knot_x = x(t0_knot_idx);

dummy = optimizableVariable('dummy', [0, 1]);
t0_optim_vars(n_t0_knots) = dummy;
for ii = 1:n_t0_knots
    name = ['t0Knot', num2str(ii)];
    t0_knot_y_ctr = t0_raw(t0_knot_idx(ii));
    range = t0_knot_y_ctr*t0_range;
    t0_optim_vars(ii) =  optimizableVariable(name, range);
end

if ~isempty(fwhm_fixed)
    self.fwhm = fwhm_fixed;
    fwhm_knot_x = [];
    fwhm_optim_vars = [];
else
    fwhm_knot_idx = get_knot_indices(p, n_fwhm_knots);
    fwhm_knot_x = x(fwhm_knot_idx);
    
    fwhm_optim_vars(n_fwhm_knots) = dummy;
    for ii = 1:n_fwhm_knots
        name = ['fwhmKnot', num2str(ii)];
        fwhm_optim_vars(ii) =  optimizableVariable(name, fwhm_range);
    end
end

% Model properties.
prop.fokin = self;
prop.x = x;
prop.t0_knot_x  = t0_knot_x;
prop.fwhm_knot_x = fwhm_knot_x;

optim_vars = [t0_optim_vars, fwhm_optim_vars];
fun = @(optim_vars_tbl)log10(do_calc_t0_fwhm(prop, optim_vars_tbl));
BO = bayesopt(fun, optim_vars, controls{:});

% Get results.
tbl = BO.XAtMinEstimatedObjective;
tbl_data = table2array(tbl);
t0_knot_y = tbl_data(1:n_t0_knots);
t0 = spline_interp(t0_knot_x, t0_knot_y, x);
self.t0 = t0;

% Calculete relative parameters.
t0_rel_params = zeros(1, n_t0_knots);
for ii = 1:n_t0_knots
    val = t0_knot_y(ii);
    range = t0_optim_vars(ii).Range;
    t0_rel_params(ii) = (val - range(1))/(range(2) - range(1));
end

if isempty(fwhm_fixed)
    fwhm_knot_y = tbl_data((n_t0_knots+1):end);
    fwhm = spline_interp(fwhm_knot_x, fwhm_knot_y, x);
    self.fwhm = fwhm;
    
    fwhm_rel_params = zeros(1, n_fwhm_knots);
    for ii = 1:n_fwhm_knots
        val = fwhm_knot_y(ii);
        range = fwhm_optim_vars(ii).Range;
        fwhm_rel_params(ii) = (val - range(1))/(range(2) - range(1));
    end
else
    fwhm = fwhm_fixed;
    fwhm_rel_params = [];
end

rel_params = [t0_rel_params, fwhm_rel_params];
end


%%
function idx = get_knot_indices(p, n_knots)
if (n_knots == 1)
    idx = round(p/2);
else
    idx = round(1:((p-1)/(n_knots-1)):p);
end
end


%%
function MSE_total = do_calc_t0_fwhm(prop, optim_vars_tbl)
fokin = prop.fokin;
x = prop.x;
t0_knot_x = prop.t0_knot_x;
fwhm_knot_x = prop.fwhm_knot_x;
n_t0_knots = numel(t0_knot_x);

knot_y = table2array(optim_vars_tbl);
t0_knot_y = knot_y(1:n_t0_knots);
fokin.t0 = spline_interp(t0_knot_x, t0_knot_y, x);

if ~isempty(fwhm_knot_x)
    fwhm_knot_y = knot_y((n_t0_knots + 1):end);
    fokin.fwhm = spline_interp(fwhm_knot_x, fwhm_knot_y, x);
end

fokin.optimize()
res = fokin.results();
MSE_total = res.MSE_total;
end


%%
function y = spline_interp(knot_x, knot_y, x)
if (numel(knot_x) == 1)
    y = knot_y;
else
    y = spline(knot_x, knot_y, x);
end
end
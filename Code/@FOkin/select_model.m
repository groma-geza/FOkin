function selected = select_model(self, lambda, omega, maxiter, info,...
    method, warm)
arguments
    self
    lambda double {mustBeNonnegative}
    omega double {mustBeNonnegative, mustBeLessThanOrEqual(omega,1)}
    maxiter double {mustBePositive, mustBeInteger}
    info(1,1) struct = struct;
    method {mustBeMember(method, {'kCV', 'RCVnv'})} = 'RCVnv';
    warm(1,1) logical = false;
end

varargin = {warm};

assert((numel(lambda) <= 2),...
    'The number of elements in argument ''lambda'' has to be 1 or 2.')
assert((numel(omega) <= 2),...
    'The number of elements in argument ''omega'' has to be 1 or 2.')
assert(~(isscalar(lambda) && isscalar(omega)),...
    ['Either argument ''lambda'' or ''omega'' has to be ',...
        'a vector of 2 elements.'])

% Create a deep copy of the options to make the following changes temporal
% and for the output.
opt = self.options.copy();
opt.bo_UseParallel = false; % Avoid nested parpools with the functions
                            % executing cross-valitation.

% Get Name, Value input parameter pairs for bayesopt() from the properties
% of opt.
controls = get_BO_controls(opt, maxiter);

% Close figures remained from previous interrupted runs.
close(findobj('Tag', 'bayesopt.ObjectiveModel'));
close(findobj('Tag', 'bayesopt.AcqFcn'));

crossval = @(lambda, omega, varargin)logCV(self, lambda, omega,...
            method, varargin{:});
        
self.select_model_running = true;
tstart = tic;
if isscalar(lambda)
    omega_optvar = optimizableVariable('omega', omega, 'Transform', 'log');
    fun = @(x)crossval(lambda, x.omega, varargin{:});
    BO = bayesopt(fun, omega_optvar, controls{:});
    selected = FOkin.SelectedOmega;
elseif isscalar(omega)    
    lambda_optvar = optimizableVariable('lambda', lambda, 'Transform', 'log');
    fun = @(x)crossval(x.lambda, omega, varargin{:});
    BO = bayesopt(fun, lambda_optvar, controls{:});
    selected = FOkin.SelectedLambda;
else
    lambda_optvar = optimizableVariable('lambda', lambda, 'Transform', 'log');
    omega_optvar = optimizableVariable('omega', omega, 'Transform', 'log');
    fun = @(x)crossval(x.lambda, x.omega, varargin{:});
    BO = bayesopt(fun, [lambda_optvar, omega_optvar], controls{:});
    selected = FOkin.SelectedLambdaOmega;
    selected.sections = FOkin.Sections(selected, []);
end
runtime = toc(tstart);
switch method
    case 'kCV'
        spec = [num2str(opt.cv_nfold), '-fold CV'];
    case 'RCVnv'
        spec = ['RCV(n_v) ', num2str(opt.cv_nrep, '%.2E'), ' repetitions'];        
end

selected.name = self.name;
selected.timestamp = char(datetime);
selected.lambda_range = lambda;
selected.omega_range = omega;
selected.BO = BO;
selected.CV_spec = spec;
selected.maxiter = maxiter;
selected.runtime = runtime;

selected.info.weight = self.weight_;
selected.info.t0 = self.t0_;
selected.info.fwhm = self.fwhm_;
selected.info.start = self.start_;

if (numel(varargin) >= 1)
    warm = varargin{1};
else
    warm = false;
end
if isempty(warm)
    warm = false;
end
selected.info.warm = warm;

selected.info.opt.fokin = obj2struct(self.options_);
selected.info.opt.elnet = obj2struct(self.optimizer_.options);
selected.info.previous = info;

if isscalar(lambda)
    [~, selected.best_omega, selected.best_value, ~,...
        selected.best_omega_1STD] = FOkin.process_selected(selected);
elseif isscalar(omega)    
    [selected.best_lambda, ~, selected.best_value,...
        selected.best_lambda_1STD] = FOkin.process_selected(selected);
else
    [selected.best_lambda, selected.best_omega, selected.best_value] =...
        FOkin.process_selected(selected);
end

% Avoid confusing new figures with a previous ones.
h_model = findobj('Tag', 'bayesopt.ObjectiveModel');
if ~isempty(h_model)
    h_model.Tag = 'bayesopt.ObjectiveModel previous';
end
h_acq = findobj('Tag', 'bayesopt.AcqFcn');
if ~isempty(h_acq)
h_acq.Tag = 'bayesopt.AcqFcn previous';
end

self.build_results('select_model')
self.select_model_running = false;

% Save a copy of of the FOkin object itself with all the data for further
% calculations.
selected.fokin = self.copy();

end


%%
function [log10_meanMSPE, constraints, user_data] = logCV(self, lambda,...
    omega, method, varargin)
[meanMSPE, constraints, user_data] = self.CV(lambda, omega, method,...
    varargin{:});
log10_meanMSPE = log10(meanMSPE);
end


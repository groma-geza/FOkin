function [meanMSPE, constraints, user_data] = CV(self, lambda, omega,...
    method, warm)
% The signature corresponds to that required by the 'bayesopt' function
% in the Statistics and Machine Learning Toolbox of MATLAB
arguments
    self
    lambda(1,1) double {mustBeNonnegative} 
    omega(1,1) double {mustBeNonnegative,...
            mustBeLessThanOrEqual(omega,1)}
    method {mustBeMember(method, {'kCV', 'RCVnv'})} = 'RCVnv';
    warm(1,1) logical = false;
end

constraints = []; % dummy outpout argument (for bayesopt)

if ~self.ready
    error('The FOkin object is not ready for optimization.')
end
self.results_ = []; % delete earliar results

options = self.options_;
optimizer = self.optimizer;
optimizer.reset() % Ensures cold restart in every step if the parfor loop
                  % runs in parallel mode or in the first step if does it
                  % in sequential mode, regardless og the value of 'warm'.

group_lasso = options.group_lasso;
if group_lasso && ~(optimizer.supports_group)
    error('The optimizer does not support group fit.')
end

optimizer.lambda = lambda;
optimizer.omega = omega;
optimizer.name = self.name_;

p = self.p;
A = self.A;
bfit = self.bfit;
multi_A = self.multi_A;

simple = ((p == 1) || group_lasso);

if multi_A
    weightfit = self.weightfit;
else
    weightfit = cell(1, p);
    for jj = 1:p
        weightfit{jj} = self.weightfit(:,jj);
    end
end

if ~multi_A
    % Reformat A and bfit according to a multi A arrangement.
    A_orig = A;
    bfit_orig = bfit;
    A = cell(1, p);
    bfit = cell(1, p);
    for jj = 1:p
        A{jj} = A_orig;
        bfit{jj} = bfit_orig(:,jj);
    end
end

len = NaN(1, p);
for jj = 1:p
    len(jj) = length(bfit{jj});
end

num_par_workers = options.num_par_workers;
if ((num_par_workers > 0) && isempty(gcp('nocreate')))
    if isinf(num_par_workers)
        parpool(options.parpool_spec);
    else
        parpool(options.parpool_spec, num_par_workers);
    end
end

switch method
    case 'kCV'
        [meanMSPE, user_data] = kCV(self, A, bfit, weightfit, len,...
            simple, warm);
    case 'RCVnv'
        [meanMSPE, user_data] = RCVnv(self, A, bfit, weightfit, len,...
            simple, warm);
end


self.runtime = [];
self.numiter = [];
self.objval = [];

self.build_results(method)
end

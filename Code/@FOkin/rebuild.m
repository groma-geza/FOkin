function rebuild(self)
self.ready = false;
self.x_norm = [];
self.objval = [];
self.runtime = [];
self.results_ = [];

params = {self.weight_, self.t0_, self.fwhm_, self.start_...
    self.options_, self.optimizer_}; % self.name_ can be empty
if any(cellfun(@isempty, params))
    return
end
p = self.p;
t0 = self.t0_;
fwhm = self.fwhm_;
start = self.start_;
weight = self.weight_;
options = self.options;

if isscalar(t0) && isscalar(fwhm) && isscalar(start)
    self.multi_A = false;
else
    if isscalar(t0) % 1xp matrix
        t0 = ones(1, p)*t0;
    end
    
    if isscalar(fwhm) % 1xp matrix
        fwhm = ones(1, p)*fwhm;
    end
    
    if isscalar(start) % 1xp matrix
        start = ones(1, p)*start;
    end
    
    % If all elements of t0, fwhm and start are identical
    % each of them can be represented by a scalar. If all columns of weight
    % are identical it can be represented by the first one. If all these
    % conditions are satisfied build a single A matrix, otherwise multi
    % ones.
    t01 = t0(1);
    cond1 = all(t0 == t01);
    fwhm1 = fwhm(1);
    cond2 = all(fwhm == fwhm1);
    start1 = start(1);
    cond3 = all(start == start1);
    cond4 = ~any(diff(weight, 1, 2), 'all');
    
    if cond1 && cond2 && cond3 && cond4
        t0 = t01;
        fwhm = fwhm1;
        start = start1;
        self.multi_A = false;
    else
        self.multi_A = true;
    end
end
self.t0_ = t0;
self.fwhm_ = fwhm;
self.start_ = start;

time = self.time;
weight = self.weight;

% Build tau
ext_lo = options.extension_lo;
ext_hi = options.extension_hi;
log_tau_first = log10(time(1)) - ext_lo;
log_tau_last = log10(time(end)) + ext_hi;

step = 1/options.n_tau_decade;
logtau = log_tau_first:step:log_tau_last;
if (logtau(end) ~= log_tau_last)
    logtau = [logtau, log_tau_last];
end
tau = 10.^logtau;

% Extend the vector of tau with inf for a constant in the measurement
% matrix
tau = [tau, inf];
self.tau = tau;
self.n = length(tau);

% Build the fitting problem
if self.multi_A
    self.A = cell(1, p);
    time2 = cell(1, p);
    self.b = cell(1, p);
    self.bfit = cell(1, p);
    self.weightfit = cell(1,p);
    
    for ii=1:p
        time2{ii} = time(start(ii):end);
        test = (time2{ii} <= 0);
        assert(~any(test), ['Delay contains negative or zero value.',...
                'Consider to increase the start parameter.'])
        
        A_ii = FOkin.design_matrix(time2{ii}, t0(ii), tau, fwhm(ii));
        w_ii = weight(start(ii):end,ii);
        self.weightfit{ii} = w_ii;
        % A weighted fit can be achieved by multiplying both self.b and the
        % colums of self.A by the square root of weightfit. This modification
        % has to be taken into account on taking the results of the fit.
        sqrt_wii = sqrt(w_ii);
        A_ii = A_ii.*sqrt_wii;
        self.A{ii} = A_ii;
        
        b_ii = self.data_norm(:,ii);
        self.b{ii} = b_ii;
        bfit_ii = b_ii(start(ii):end).*sqrt_wii;
        self.bfit{ii} = bfit_ii;    
    end
else
    time2 = time(start:end);
    test = (time2 <= 0);
    assert(~any(test), ['Delay contains negative or zero value.',...
            'Consider to increase the start parameter.'])
    
    A = FOkin.design_matrix(time2, t0, tau, fwhm);
    weightfit = weight(start:end,:);
    self.weightfit = weightfit;
    % A weighted fit can be achieved by multiplying both self.b and the
    % colums of self.A by the square root of weightfit. This modification
    % has to be taken into account on taking the results of the fit.
    sqrt_weightfit = sqrt(weightfit);
    % All the columns of sqrt_weightfit are identical, use the first one
    % MATLAB earlier than 2016b needs bsxfun here.
    A = A.*sqrt_weightfit(:,1);    
    self.A = A;
    
    b = self.data_norm;
    self.b = b;
    bfit = b(start:end,:).*sqrt_weightfit;
    self.bfit = bfit;
end

self.select_model_running = false;
self.ready = true;
end



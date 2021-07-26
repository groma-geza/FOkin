function build_results(self, operation)
% Collect results of optimization.    

switch operation
    case {'select_model', 'no operation'}
        level = 0;
    case {'kCV', 'RCVnv'}
        level = 1;
    case 'optimize'
        level = 2;
    otherwise
        error('Invalid input argument.')
end

multi_A = self.multi_A;
n = self.n -1;
p = self.p;
norm = self.norm;
grouped = self.options_.group_lasso;

res = FOkin.FOkinResults;
res.level = level;
res.norm_orig = norm;
res.optimizer_infotext1 = self.optimizer_.infotext1();
res.optimizer_infotext2 = self.optimizer_.infotext2();
res.optimizer_class = class(self.optimizer_);

res.name = self.name;
res.options = self.options_.copy();
res.operation = operation;
res.timestamp = char(datetime);
res.time = self.time;
res.group_param = self.group_param;
res.tau = self.tau(1:(end-1));
res.norm = norm;

if (level > 0)
    x_norm = self.x_norm;
    if ~multi_A
        fit_total = self.A*x_norm;
    end
    res.x0 = zeros(1, p);
    res.x = zeros(n, p);
end

for ii = p:-1:1 %backwards for preallocation
    if (level > 0)
        x = x_norm(:,ii);
    end
    if multi_A
        b = self.b{ii};
        bfit = self.bfit{ii};
        if (level > 0)
            fit = self.A{ii}*x;
        end
        start = self.start(ii);
        t0 = self.t0(ii);
        fwhm = self.fwhm(ii);
        weightfit = self.weightfit{ii};               
    else
        b = self.b(:,ii);
        bfit = self.bfit(:,ii);
        if (level > 0)
            fit = fit_total(:,ii);
        end
        start = self.start;
        t0 = self.t0;
        fwhm = self.fwhm;
        weightfit = self.weightfit(:,ii);
    end
    
    % Compensate the effect of the method applied for acheaving a weighted
    % fit in self.rebuild()
    sqrt_weightfit = sqrt(weightfit);
    bfit = bfit./sqrt_weightfit;
    if (level > 0)
        fit = fit./sqrt_weightfit;
        
        residual = bfit - fit;
        residual_total(:,ii) = residual;
        weight_total(:,ii) = weightfit;
        MSE = wmse(residual, weightfit);
        
        res.x(:,ii) = x(1:(end-1));
        res.x0(ii) = x(end);
    end
    fd = FOkin.Fitdata;
    fd.group_param = self.group_param(ii);
    fd.b = b;
    fd.t0 = t0;
    fd.fwhm = fwhm;
    fd.start = start;    
    fd.bfit = bfit;
    fd.weightfit = weightfit;
    if (level > 0)
        fd.fit = fit;
        fd.residual = residual;
    end
    % self.objval is obtained from the fit of the original
    % normalized and weighted data.
    if (~grouped && (level == 2))
        fd.objval = self.objval(ii);
    end
    if (level > 0)
        fd.MSE = MSE;
    end
    fitdata(ii) = fd;    
end

res.fitdata = fitdata;

if (level > 0)
    res.lambda = self.optimizer_.lambda;
    res.omega = self.optimizer_.omega;
end

if (level == 2)
    res.numiter = self.numiter;
    res.runtime = self.runtime;
    if ~isempty(self.runtime)
        res.runtime_total = sum(self.runtime);
    end
end

% self.objval is obtained from the fit of the original
% normalized and weighted data.
if (grouped && (level == 2))
    res.group_objval = self.objval;
end

if (level > 0)
    MSE_total = wmse(residual_total, weight_total);
    res.MSE_total = MSE_total;
end

res.info.name = self.name;
res.info.weight = self.weight_;
res.info.t0 = self.t0_;
res.info.fwhm = self.fwhm_;
res.info.start = self.start_;
res.info.opt.fokin = obj2struct(self.options_);
res.info.opt.elnet = obj2struct(self.optimizer_.options);

self.results_ = res;
end

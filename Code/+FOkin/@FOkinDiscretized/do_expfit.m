function do_expfit(self)
% Fits exponentials to the data using the discretized tau values which are
% not neglected.

time = self.time;
group_param = self.group_param;
p = length(group_param);
fitdata = self.fitdata;
% Collegth length of data in the elements of fitdata
m = NaN(1, p);
for ii = 1:p
    m(ii) = length(fitdata(ii).bfit);
end

tau_in = self.discrete.average_tau';

scale = tau_in;

prm0 = ones(size(scale));
fitfun = @(prm)expfitfun(prm, scale, time, fitdata, m);
[prm, ~, ~, exitflag, output] = lsqnonlin(fitfun, prm0);

if exitflag >= 0
    [residual_total, DADS, fit, MSE, MSE_total] =...
        expfitfun(prm, scale, time, fitdata, m);
    
    tau = (prm.*scale)';
    abs_amplitude = max(abs(DADS), [], 2);
    abs_max_amplitude = max(abs_amplitude);
    rel_amplitude = abs_amplitude/abs_max_amplitude;

    % Sort rel_amplitude, DADS and tau in descending rel_amplitude.
    [rel_amplitude, idx] = sort(rel_amplitude,'descend' );
    DADS = DADS(idx,:);
    tau = tau(idx);

    expfit = FOkin.Expfit;
    expfit.tau = tau;
    expfit.DADS = DADS;
    expfit.rel_amplitude = rel_amplitude;
    expfit.abs_max_amplitude = abs_max_amplitude;
    expfit.fit = fit;
    expfit.residual = cell(1, p);
    ptr = 1;
    for ii = 1:p
        new_ptr = ptr + m(ii);
        expfit.residual{ii} = residual_total(ptr:(new_ptr-1));
        ptr = new_ptr;
    end
    expfit.MSE = MSE;
    expfit.MSE_total = MSE_total;    
    self.discrete.expfit = expfit;
else
    error(output.message)
end
end

%%
function [residual_total, DADS, fit, MSE, MSE_total] = ...
    expfitfun(prm, scale, t, fitdata, m)

tau = prm.*scale;
n = length(tau);
p = length(fitdata);
DADS = NaN(n, p);
fit = cell(1, p);
MSE = NaN(1, p);

m_total = sum(m);
residual_total = NaN(1, m_total);
weight_total = NaN(1, m_total);
ptr = 1;
for ii = 1:p
    % Solve the linear optimization problem for the amplitudes
    % Numerical Recipes in C (1989) pp. 528-534
    fd_ii = fitdata(ii);
    t0 = fd_ii.t0;
    fwhm = fd_ii.fwhm;
    start = fd_ii.start;
    % The columns of B are the base functions.    
    B = FOkin.design_matrix(t(start:end), t0, tau, fwhm);
    bfit = fd_ii.bfit;
    weightfit = fd_ii.weightfit;
    DADS(:,ii) = (B'*B)\(B'*bfit); % vector of amplitudes
    fit_ii = B*DADS(:,ii);
    fit{ii} = fit_ii;
    residual = (bfit - fit_ii).*sqrt(weightfit);
    new_ptr = ptr + m(ii);
    residual_total(ptr:(new_ptr-1)) = residual; % weighted residual
    weight_total(ptr:(new_ptr-1)) = weightfit;
    MSE(ii) = sum(residual.^2)/sum(weightfit);
    ptr = new_ptr;
end
MSE_total = sum(residual_total.^2)/sum(weight_total);
end
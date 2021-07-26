function [best_lambda, best_omega ,best_value, best_lambda_1STD,...
    best_omega_1STD] = process_selected(selected, lambda_slice_width,...
    omega_slice_width, show)
arguments
    selected(1,1) FOkin.SelectedModel
    lambda_slice_width(1,1) double {mustBePositive} = Inf;
    omega_slice_width(1,1) double {mustBePositive} = Inf;
    show(1,1) {mustBeMember(show, [0, 1, 2])} = 0;
end

BO = selected.BO;
name = selected.name;
spec = selected.CV_spec;

min_data_points = 20;
predict_size = 1000;
    
lambda_only = isscalar(selected.omega_range);
omega_only = isscalar(selected.lambda_range);

if lambda_only
    lim_lambda_lo = BO.VariableDescriptions(1, 1).Range(1);
    lim_lambda_hi = BO.VariableDescriptions(1, 1).Range(2);
    lim_omega_lo = -Inf;
    lim_omega_hi = Inf;
    lambda_trace = BO.XTrace.lambda;
    omega_single = selected.omega_range;
    omega_trace = omega_single*ones(size(lambda_trace));
    omega_slice_width = Inf;
elseif omega_only
    lim_lambda_lo = -Inf;
    lim_lambda_hi = Inf;
    lim_omega_lo = BO.VariableDescriptions(1, 1).Range(1);
    lim_omega_hi = BO.VariableDescriptions(1, 1).Range(2);
    omega_trace = BO.XTrace.omega;   
    lambda_single = selected.lambda_range;
    lambda_trace = lambda_single*ones(size(omega_trace));
    lambda_slice_width = Inf;
else
    lim_lambda_lo = BO.VariableDescriptions(1, 1).Range(1);
    lim_lambda_hi = BO.VariableDescriptions(1, 1).Range(2);
    lim_omega_lo = BO.VariableDescriptions(1, 2).Range(1);
    lim_omega_hi = BO.VariableDescriptions(1, 2).Range(2);
    lambda_trace = BO.XTrace.lambda;
    omega_trace = BO.XTrace.omega;
end

obj_trace = BO.ObjectiveTrace;

if ~omega_only
    log_first = log10(lim_lambda_lo);
    log_last = log10(lim_lambda_hi);
    step = (log_last - log_first)/(predict_size - 1);
    loglambda = log_first:step:log_last;
    lambda_vals = (10.^loglambda)';
else
    lambda_vals = lambda_single*ones(predict_size, 1);
end
if ~lambda_only
    log_first = log10(lim_omega_lo);
    log_last = log10(lim_omega_hi);
    step = (log_last - log_first)/(predict_size - 1);
    logomega = log_first:step:log_last;
    omega_vals = (10.^logomega)';
else
    omega_vals = omega_single*ones(predict_size, 1);
end

[best_point, best_log_value] = BO.bestPoint('Criterion', 'min-mean');
% The return value of BO.bestPoint is stochastic, use it only once for
% an object.
if ~omega_only
    best_lambda = selected.best_lambda;
    if isempty(best_lambda)
        best_lambda = best_point.lambda;
    end
else
    best_lambda = lambda_single;
end
if ~lambda_only
    best_omega = selected.best_omega;
    if isempty(best_omega)
        best_omega = best_point.omega;
    end
else
    best_omega = omega_single;
end

best_value = selected.best_value;
if isempty(best_value)
    best_value = 10^best_log_value;
end

if ((show == 2) && ~lambda_only && ~omega_only)
    BO.plot(@plotObjectiveModel);
    hfig = gcf;
    hfig.Tag = 'bayesopt.ObjectiveModel previous';
    hax = gca;
    hax.XLabel.String = '\lambda';
    hax.YLabel.String = '\omega';
    hax.ZLabel.String = 'meanMSE';
    hlin = findobj(hax, 'DisplayName', 'Next point');
    delete(hlin)
    hlin = findobj(hax, 'DisplayName', 'Model minimum feasible');
    delete(hlin)
      
    hlin1 = findobj(hax, 'DisplayName', 'Observed points');
    hlin1.ZData = 10.^hlin1.ZData;
    hlin1.MarkerSize = 5;
    hsurf = findobj(hax, 'DisplayName', 'Model mean');
    xdat = hsurf.XData;
    ydat = hsurf.YData;
    zdat = 10.^hsurf.ZData;
    delete(hsurf)
    hcont = findobj(hax, 'Type', 'Contour');
    delete(hcont)
    hax.ZScale = 'log';
    hold on
    hsurf = surfc(hax, xdat, ydat, zdat, 'FaceColor', 'r', 'FaceAlpha',...
        .5, 'LineStyle', '-');
    hlin2 = line([best_lambda, best_lambda], [best_omega, best_omega],...
        zlim, 'Color','c', 'LineStyle', ':',  'LineWidth', 3);

    txt = {[name, '   Objective function model   ', spec],...
        ['\lambda_{@min} = ', num2str(best_lambda, '%.2E'), '   ',...
        '\omega_{@min} = ', num2str(best_omega, '%.2E'),  '   ',...
        'value_{@min} = ', num2str(best_value, '%.4E')],  ''};
      
    title(txt, 'FontWeight', 'normal')
    
    legend([hlin1, hsurf(1), hlin2], 'Observed points', 'Model mean',...
        'Position of minimum')
    
    BO.plot(@plotAcquisitionFunction)
    hfig = gcf;
    hfig.Tag = 'bayesopt.AcqFcn previous';    
    hax = gca;
    hax.XLabel.String = '\lambda';
    hax.YLabel.String = '\omega';
    hax.ZLabel.String = 'Expected improvement';
    txt = {[name, '   Expected improvement', '   ', spec], ''};
    title(txt, 'FontWeight', 'normal')
end

%%
% BEGIN nested function
function best_p1_1STD = get_minimum_single(tbl, p1, p1_trace, best_p1,...
        p1_label, p2_slice_width, lim_p2_lo, lim_p2_hi, p2_trace, best_p2)

[objective, FSD, noise_sigma] = BO_mean_err(BO, tbl);

is_lambda = strcmp(p1_label, 'lambda');
if is_lambda
    p2_label = 'omega';
else
    p2_label = 'lambda';
end

if isinf(p2_slice_width)
    p1_mask = p1_trace;
    obj_mask = obj_trace;
else
    % best_p2 can not be 0
    p2_lo = 10^(log10(best_p2) - p2_slice_width);
    p2_hi = 10^(log10(best_p2) + p2_slice_width);
    if p2_lo < lim_p2_lo
        p2_lo = lim_p2_lo;
        p2_hi = 10^(log10(p2_lo) + 2*p2_slice_width);
    end
    if p2_hi > lim_p2_hi
        p2_hi = lim_p2_hi;
        p2_lo = 10^(log10(p2_hi) - 2*p2_slice_width);
    end
    
    % Select the traces around best_p2
    mask = ((p2_trace >= p2_lo) & (p2_trace <= p2_hi));
    mask_size = sum(mask);
    if (~isinf(p2_slice_width) &&(mask_size < min_data_points))
        error(['Only ', num2str(mask_size),...
            ' data point(s) fall(s) into the selected data range. ',...
            'Increase ''', p2_label, '_slice_width''!'])
    end
    p1_mask = p1_trace(mask);
    obj_mask = obj_trace(mask);
end

% Find nearest neighbours of best_p1 in p1_vals
lo = find(p1 <= best_p1);
lo = lo(end);
hi = find(p1 >= best_p1);
hi = hi(1);
if lo == hi
    if hi < predict_size
        hi = hi + 1;
    else
        lo = lo - 1;
    end
end
level_min = interp1(p1, objective, best_p1);
x = [p1(lo), p1(hi)];
v = [objective(lo) + noise_sigma, objective(hi) + noise_sigma];
level_1STD = interp1(x, v, best_p1);

% Apply 1 STD rule with sigma
idx = find(objective <= level_1STD);
invalid = false;
if ~isempty(idx)
    % Take 1 STD toward the simpler model
    if is_lambda
        idx = idx(end); % the higher is the simpler
    else
        idx = idx(1) - 1; % the lower is the simpler
    end
    if idx == predict_size
        best_p1_1STD = p1(idx);
    elseif idx == 0
        best_p1_1STD = p1(1);
    else
        % Make interpolation
        x1 = p1(idx);
        x2 = p1(idx + 1);
        y1 = objective(idx);
        y2 = objective(idx + 1);
        best_p1_1STD = x1 + (x2 - x1)*(level_1STD -y1)/(y2 -y1);
    end
    
    if is_lambda && (best_p1_1STD <= best_p1)
        invalid = true;
    end
    if ~is_lambda && (best_p1_1STD >= best_p1)
        invalid = true;
    end
else
    invalid = true;
end
if invalid
    best_p1_1STD = NaN;
end

if (show == 1)
    figure('Position', [60, 80, 800, 576])
    hl1 = loglog(p1_mask, 10.^obj_mask, 'bo', 'MarkerSize', 5,...
        'MarkerFaceColor', 'b');
    hold on
    h2 = loglog(p1, 10.^objective, 'r');
    h3 = loglog(p1, 10.^(objective + FSD), 'b');
    loglog(p1, 10.^(objective - FSD), 'b');
    h4 = loglog(p1, 10.^(objective + noise_sigma), 'c');
    loglog(p1, 10.^(objective - noise_sigma), 'c');

    xlim([p1(1)*.9, p1(end)*1.1])
    y_lim = ylim;
    h5 = line([best_p1, best_p1], y_lim, 'Color','r', 'LineStyle', ':',...
        'LineWidth', 1.5);
    h6 = line([best_p1_1STD, best_p1_1STD], y_lim, 'Color','g',...
        'LineStyle', ':',  'LineWidth', 1.5);
    line(xlim, 10.^[level_min, level_min], 'Color','r',...
        'LineStyle', ':',  'LineWidth', 1.5);
    line(xlim, 10.^[level_1STD, level_1STD], 'Color','g',...
        'LineStyle', ':',  'LineWidth', 1.5);
    hold off
    ylim(y_lim)
    xlabel(['\', p1_label])
    ylabel('mean MSE')
    
    txt = {[name, '   Objective function model   ', spec],...
        ['@\',p2_label, ' = ', num2str(best_p2,...
        '%.2E'), '   \',p1_label,...
        '_{@min} = ', num2str(best_p1, '%.2E'),...
        '   value_{@min} = ', num2str(best_value, '%.4E'),...
        '   ', '\', p1_label,...
        '_{@min + 1 STD} = ', num2str(best_p1_1STD, '%.2E')], ''};
      
    title(txt, 'FontWeight', 'normal')
    
    legend([hl1, h2, h3, h4, h5, h6], 'Observed points', 'Model mean',...
        'Model error bars', 'Noise error bars', 'Position of minimum',...
        'Position of minimum + 1 STD')
       
    % Expected improvement (plus) according to the undocumented
    % bayesoptim.expectedImprovement function of the Statistics and 
    % Machine Learning Toolbox of MATLAB.    
    GammaX = (min(objective) - objective)./FSD;
    PI = normcdf(GammaX, 0, 1);
    EI = FSD.*(GammaX.*PI + normpdf(GammaX, 0, 1));
    
    figure, semilogx(p1, EI)
    xlabel(['\', p1_label])
    ylabel('Expected improvement')
    txt = {[name, '   Expected improvement', '   ', spec], ''};
    title(txt, 'FontWeight', 'normal')
end

end
% END nested function

%%
if ~omega_only
    % Calulate minimum along lambda
    lambda = lambda_vals;
    omega = ones(predict_size, 1)*best_omega;
    tbl = table(lambda, omega);
    best_lambda_1STD = get_minimum_single(tbl,...
        lambda, lambda_trace, best_lambda, 'lambda',...
        omega_slice_width, lim_omega_lo, lim_omega_hi, omega_trace,...
        best_omega);
end

if ~lambda_only
    % Calulate minimum along omega
    omega = omega_vals;
    lambda = ones(predict_size, 1)*best_lambda;
    tbl = table(lambda, omega);
    best_omega_1STD = get_minimum_single(tbl,...
        omega, omega_trace, best_omega, 'omega',...
        lambda_slice_width, lim_lambda_lo, lim_lambda_hi, lambda_trace,...
        best_lambda);
else
    best_omega = [];
    best_omega_1STD = [];
end

if omega_only
    best_lambda = [];
    best_lambda_1STD = [];
end
end

%%
function [objective, FSD, noise_sigma] = BO_mean_err(BO, tbl)
% Returns the the model mean, model error and noise error of a BO object
% at the values of optimzation variables defined in a table

    % input - BO : BayesianOptimization object
    %       - tbl: table of lambda and omega optimzation variables
    
    % output - objective    : model mean
    %        - FSD          : model error
    %        - noise_sigma  : noise error


    if ~isempty(BO.ObjectiveFcnGP)
        noise_sigma = BO.ObjectiveFcnGP.Sigma; % works on MATLAB 2019b
    else
        noise_sigma =  BO.ObjectiveFcnModel.Sigma; % works on MATLAB 2020a
                                                   % and later
    end
[objective, sing_sigma] = BO.predictObjective(tbl);

% FSD is calculated according to the undocumented 
% bayesoptim.probabilityOfImprovement function in the Statistics and 
% Machine Learning Toolbox of MATLAB.

FSD = sqrt(max(0, sing_sigma.^2 - noise_sigma.^2));
end

function [t0, fwhm, MSE_total] = risefit(data, time, group_param,...
    opt, fokin_opt, elnetADMM_opt, elnetPDCO_opt)

% Input
% data          - mxp matrix storing the kinetic data in arbitrary units.
%                 A row of the matrix corresponds to a given group_param
%                 value and different time values and a column corresponds
%                 to a  given time value and different group_param values.
% time          - mx1 vector, containing the time
%                 points in increasing order but in arbitrary steps. An
%                 arrangement close to logarithmically equidistant point is
%                 suggested.
% group_param   - 1xp vector, containing a group
%                 parameter (typically wavelength) in increasing order but
%                 in arbitrary steps.
% opt           - FOkin.RisefitOptions object. 
% fokin_opt     - FOkinOptions object.
% elnetADMM_opt - MultiElnetADMMOptions object.
% elnetPDCO_opt - MultiElnetPDCOOptions object.

% Output t0     - 1xp vector, containing the data required by the t0
%                 property of a FOkin object.
% fwhm          - 1xp vector, containing the data required by the fwhm
%                 property of a FOkin object.
% MSE_total     - scalar, the weighted average of the total residual of
%                 fitting obtained by the returned value of t0 and fwhm.

arguments
    data double
    time double {mustBePositive}
    group_param double
    opt FOkin.RisefitOptions = FOkin.RisefitOptions;
    fokin_opt FOkinOptions = FOkinOptions;
    elnetADMM_opt MultiElnetADMMOptions = MultiElnetADMMOptions;
    elnetPDCO_opt MultiElnetPDCOOptions = MultiElnetPDCOOptions;
end

% Create a deep copy of fok_opt for changing values in it.
fokin_opt = fokin_opt.copy();
show_details = opt.show_details;
individual = opt.individual;
repeat = opt.repeat;
maxiter = opt.maxiter;
title_text = opt.title_text;
fwhm_fixed = fokin_opt.fwhm_fixed;

if show_details == 0
    fokin_opt.bo_Verbose = 0;
end
if show_details <= 2
    fokin_opt.bo_PlotFcn = [];
end

if fokin_opt.group_lasso
    % ADMM is required
    elnet = MultiElnetADMM();
    elnet.options = elnetADMM_opt;
else
    % PDCO is faster
    elnet = MultiElnetPDCO();
    elnet.options = elnetPDCO_opt;
end
elnet.lambda = opt.elnet_lambda;
elnet.omega = opt.elnet_omega;

p = numel(group_param);
t0 = NaN(1, p);
fwhm = NaN(1, p);

fok_opt_individ = fokin_opt.copy(); % deep copy
fok_opt_individ.n_t0_knots = 1;
fok_opt_individ.n_fwhm_knots = 1;

% Precedence of options
if (~isempty(fwhm_fixed) || ~individual)
    repeat = false;
end

if repeat
    steps = 2;
    fok_opt_individ.fwhm_fixed = [];
    if (numel(maxiter) == 1)
        maxiter = [maxiter, maxiter];
    end
else
    steps = 1;
end

for ii = 1:steps
    if (individual || (p == 1))
        if (ii == 1)
            rel_params = NaN(p, 2);
        else
            rel_params = NaN(p, 1);
        end
        
        for jj = 1:p
            fokin = FOkin(data(:,jj), time, group_param(jj));
            fokin.name = opt.title_text;
            elnet_individ = MultiElnetPDCO();
            elnet_individ.lambda = opt.elnet_lambda;
            elnet_individ.omega = opt.elnet_omega;
            fokin.optimizer = elnet_individ;
            
            if (~isempty(fwhm_fixed) && ~isscalar(fwhm_fixed))
                fok_opt_individ.fwhm_fixed = fwhm_fixed(jj);
            end
            fokin.options = fok_opt_individ;
            
            [t0(jj), fwhm(jj), rel_params(jj,:)] =...
                fokin.calc_t0_fwhm(maxiter(ii));
            
            if show_details
                fprintf('\n')
                disp('Relative parameters:')
                disp(rel_params(jj,:))
            end
        end
    else
        if opt.separate_norm
            maxvals = max(data, [], 1);
            norm_data = data./maxvals;
            fokin = FOkin(norm_data, time, group_param);
            fokin.name = opt.title_text;
        else
            fokin = FOkin(data, time, group_param);
            fokin.name = opt.title_text;
        end
        fokin.optimizer = elnet;
        fokin.options = fokin_opt;
        
        [t0, fwhm, rel_params] = fokin.calc_t0_fwhm(maxiter(ii));
        
        if show_details
            fprintf('\n')
            disp('Relative parameters:')
            disp(rel_params)
        end
    end
    if ((ii == 1) && (steps == 2) && (show_details >= 2))
        plot_t0_fwhm(group_param, t0, fwhm, [], rel_params,...
            [title_text, ' step1'], fokin_opt)
    end
    if repeat
        if (ii == 1)
            fok_opt_individ.fwhm_fixed = mean(fwhm);
        else
            fwhm = mean(fwhm); % ensure to store in scalar form
        end
    end
end

% Rebuild fokin for the original problem and optimize with that.
fokin = FOkin(data, time, group_param);
fokin.name = opt.title_text;
fokin.optimizer = elnet;
fokin.options = fokin_opt;
fokin.t0 = t0;
fokin.fwhm = fwhm;

fokin.optimize();
res = fokin.results();
MSE_total = res.MSE_total;
if show_details >= 2
    if (ii == 2)
        title_text = [title_text, ' step2'];
    end
    plot_t0_fwhm(group_param, t0, fwhm, MSE_total,...
        rel_params, title_text,  fokin_opt)
end
if show_details >= 3
    res.show
end
end


%%
function plot_t0_fwhm(lambda, t0, fwhm, MSE_total, rel_params,...
    title_text, fok_opt)

if ~isscalar(rel_params(:,1))
    rel_params = mean(rel_params);
    state = ' (average)';
else
    state = '';
end
        
figure('Position', [10, 380, 560, 420]);
yyaxis left
plot(lambda, t0, '.-')
ylabel(['t0 (', fok_opt.time_unit, ')'])

yyaxis right
if isscalar(fwhm)
    plot(lambda, ones(size(lambda))*fwhm, '.-');
    ylim(fok_opt.fwhm_range)
else
    plot(lambda, fwhm, '.-')
end
ylabel(['fwhm (', fok_opt.time_unit, ')'])
xlabel(fok_opt.group_param_label)

if ~isempty(MSE_total)
    title_text = [title_text, '       MSE_{total} = ',...
        num2str(MSE_total, '%.4E')];
end

title({title_text, ['Relative parameters', state, ': ',...
    num2str(rel_params)]}, 'FontWeight', 'normal');

end


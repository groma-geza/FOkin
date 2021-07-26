function  showexp(self, selection)
arguments
    self
    selection double = [];
end

if (self.norm == 1)
    normalized = false;
else
    normalized = true;
end

time = self.time;
group_param = self.group_param;
p = length(group_param);
fitdata = self.fitdata;
expfit = self.discrete.expfit;
if isempty(expfit.tau)
    return
end
tau = expfit.tau;
nexp = length(tau);
opt = self.options;
signal_label = opt.signal_label;
time_unit = opt.time_unit;
title_text1 = [self.name, '   Result of expfit   nexp = ', num2str(nexp)];
x_lim = [self.tau(1)/1.5, self.tau(end)*1.5];
for ii = 1:p
    if (isempty(selection)...
            || (~isempty(find(group_param(ii) == selection, 1))))
        figure('Position', [40, 80, 800, 376.8], 'Tag', 'expfit');
        h1 = subplot(2, 1, 1);
        h2 = subplot(2, 1, 2);
        
        fd_ii = fitdata(ii);
        b = fd_ii.b;
        start = fd_ii.start;
        time2 = time(start:end);
        fit = expfit.fit{ii};
        semilogx(h1, time, b, 'b.', time2, fit, 'r',...
            'LineWidth', 1, 'MarkerSize', 8)
        
        h1.XLim = x_lim;
        h1.XLabel.String = ['Time (', time_unit, ')'];
        txt = {signal_label};        
        if normalized
            txt{1, 2} = '(normalized)';
        end
        h1.YLabel.String = txt;
        h1.XGrid = 'on';
        h1.YGrid = 'on';
        h1.Title.String = {title_text1, ['   wavelength = ',...
            num2str(group_param(ii)), ' nm   MSE = ',...
            num2str(expfit.MSE(ii), '%.4E'),...
            '   MSE_{total} = ', num2str(expfit.MSE_total, '%.4E')]};
        h1.Title.FontWeight = 'normal';
        h1.Title.FontSize = 9;
        
        semilogx(h2, time2, expfit.residual{ii}, '.-r', 'LineWidth', 1)
        h2.XLim = x_lim;
        h2.XLabel.String = ['Time (', time_unit, ')']; 
        txt = {'Residual'};
        if normalized
            txt{1, 2} = '(normalized)';
        end
        h2.YLabel.String = txt;
        h2.XGrid = 'on';
        h2.YGrid = 'on';
        drawnow
    end
end

figure('Position', [860, 80, 560, 420])
h3 = axes;
h3.NextPlot = 'add';
h3.ColorOrder = opt.thread_colors;
h3.LineStyleOrder = {'-','--',':'};

s_param = (1:500)/500*(max(group_param) - min(group_param)) +...
    min(group_param);

for ii=1:nexp
    DADS_ii =  expfit.DADS(ii,:);
    DADS_spline = csaps(group_param, DADS_ii, opt.smooth_val);
    s_DADS = fnval(DADS_spline, s_param);
    coi = h3.ColorOrderIndex;
    loi = h3.LineStyleOrderIndex;
    if ~isinf(tau(ii))
        plot(h3, group_param, DADS_ii, '.', 'MarkerSize', 8)
        h3.ColorOrderIndex = coi;
        h3.LineStyleOrderIndex = loi;
        plot(h3, s_param, s_DADS, 'LineWidth', .5)
    else
        plot(h3, group_param, DADS_ii, 'k.', 'MarkerSize', 8)
        plot(h3, s_param, s_DADS, 'k', 'LineWidth', .5)
        h3.ColorOrderIndex = coi;
        h3.LineStyleOrderIndex = loi;        
    end
end

line(xlim, [0, 0], 'Color', 'k', 'LineStyle', ':')
xlabel(opt.group_param_label);
txt = 'Amplitude';
if normalized
    txt = [txt, ' (nomalized)'];
end
ylabel(txt)
title([title_text1, '   MSE_{total} = ',...
    num2str(expfit.MSE_total, '%.4E')], 'FontWeight','normal',...
    'FontSize', 9)

self.show_table('expfit')
end
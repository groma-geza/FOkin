function show(self, selection)
% selection  - an array of values of self.group_param for which display is
%               required.
arguments
    self
    selection double = [];
end

assert((self.level > 0), 'No results to show.')

if (self.norm == 1)
    normalized = false;
else
    normalized = true;
end

x_total = self.x;
x0 = self.x0;
time = self.time;
group_param = self.group_param;
p = length(group_param);
tau = self.tau;
opt = self.options;
signal_label = opt.signal_label;
time_unit = opt.time_unit;
group_param_name = opt.group_param_name;

figure
h = axes;
semilogx(tau, x_total, tau(end)*1.2, x0, 'k.')
x_lim = [tau(1)/1.5, tau(end)*1.5];
h.XLim = x_lim;
min_x_total = min(min(x_total));
max_x_total = max(max(x_total));
range = max_x_total - min_x_total;
if range ~= 0
    h.YLim = [min_x_total - range/20, max_x_total + range/20];
end

h.XLabel.String = ['\tau (', time_unit, ')'];
txt = {'Amplitude'};
if normalized
    txt{1, 2} = '(normalized)';
end
h.YLabel.String = txt;
txt = [self.name, '   Result of fit'];
if p > 1
    txt = [txt, ' at different ' , group_param_name, ' values'];
end
txt = {txt, self.short_header(), ''};
h.Title.String = txt;
h.Title.FontSize = 9;
h.Title.FontWeight = 'normal';
grid on

if ((p > 1) && (~all(x_total == 0, 'all')))
    figure
    hax = self.plottau2D(x_total);
    txt = [self.name, '   Result of fit'];    
    if normalized
        txt = [txt, '   ', '(normalized)'];
    end
    txt = {txt, self.short_header(), ''};
    hax.Title.String = txt;
    hax.Title.FontSize = 9;
    hax.Title.FontWeight = 'normal';
end

for ii = 1:p
    x = x_total(:,ii);
    fd = self.fitdata(ii);
    b = fd.b;
    fit = fd.fit;
    residual = fd.residual;
    start = fd.start;    
    time2 = time(start:end);
    
    if (isempty(selection)...
            || (~isempty(find(group_param(ii) == selection, 1))))
        figure('Position', [850, 40, 800, 650]);
        h1 = subplot(3, 1, 1);
        h2 = subplot(3, 1, 2);
        h3 = subplot(3, 1, 3);
        
        %         semilogx(h1, tau, all_x, 'LineWidth', 1);
        hbr = bar(h1, tau, x, 'BarWidth', 0, 'LineWidth', 1.5,...
            'EdgeColor', 'blue');
        hbr.BaseLine.Color = 'none';
        
        h1.NextPlot = 'add';
        h1.XGrid = 'on';
        h1.YGrid = 'on';
        
        h1.XScale = 'log';
        h1.XLim = x_lim;
        min_x = min(x);
        max_x = max(x);
        range = max_x - min_x;
        
        hbr_const = bar(h1, tau(end)*1.2, x0(ii),  'BarWidth', 0,...
            'LineWidth', 1.5, 'EdgeColor', 'k');
        hbr_const.BaseLine.Color = 'none';
        
        if range ~= 0
            h1.YLim = [min_x - range/20, max_x + range/20];
            hbr.BaseValue = range/100;
        end
        txt = [self.long_header(ii), {''}];
        h1.Title.String = txt;
        h1.Title.FontWeight = 'normal';
        h1.Title.FontSize = 9;
        h1.XLabel.String = ['\tau (', time_unit, ')'];
        txt = {'Amplitude'};
        if normalized
            txt{1, 2} = '(normalized)';
        end
        h1.YLabel.String = txt;
        h1.NextPlot = 'add';
        h1.YLimMode = 'manual';
        
        semilogx(h2, time, b, 'b.', time2, fit, 'r',...
            'LineWidth', 1, 'MarkerSize', 8)
        
        h2.XLim = x_lim;
        h2.XLabel.String = ['Time (', time_unit, ')'];
        txt = {signal_label};        
        if normalized
            txt{1, 2} = '(normalized)';
        end
        h2.YLabel.String = txt;
        h2.XGrid = 'on';
        h2.YGrid = 'on';
        %         legend({'experimental', 'fit'}, 'EdgeColor', 'white',...
        %               'Location', 'northwest')
        
        semilogx(h3, time2, residual, '.-r', 'LineWidth', 1)
        h3.XLim = x_lim;
        h3.XLabel.String = ['Time (', time_unit, ')'];
        txt = {'Residual'};
        if normalized
            txt{1, 2} = '(normalized)';
        end
        h3.YLabel.String = txt;
        h3.XGrid = 'on';
        h3.YGrid = 'on';
        drawnow
    end
end
end
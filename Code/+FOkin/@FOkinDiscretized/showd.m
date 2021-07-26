function showd(self)

if (self.norm == 1)
    normalized = false;
else
    normalized = true;
end

x = self.x;
if all(all(x == 0))
    return
end

disc = self.discrete;
param = self.group_param;
opt = self.options; 
threads = disc.threads;
k = length(threads);

if ((length(param) > 1) && (k~=0))
    % Display results
    figure, hax = self.plottau2D(x);
    figure('Position', [660, 230, 1000, 510]);
    co = opt.thread_colors;
    lso = {'-','--',':'};    
    h2 = subplot('Position', [0.06, 0.1, 0.42, 0.82]);
    h3 = subplot('Position', [0.56, 0.1, 0.42, 0.82]);

    hax.NextPlot = 'add';
    h2.YAxis.Scale = 'log';
    h2.NextPlot = 'add';    
    h3.NextPlot = 'add';
    h2.ColorOrder = co;
    h3.ColorOrder = co;
    h2.LineStyleOrder = lso;
    h3.LineStyleOrder = lso;
    
    idx_x0 = [];
    for ii=1:k
        feature = threads(ii);       
        f_param = feature.param;
        f_tau = feature.tau;
        if isinf(f_tau)
            % thread storing x0
            idx_x0 = ii;
        else
            f_val = feature.val;
            s_param = (1:500)/500*(max(f_param) - min(f_param)) +...
                min(f_param);
            tau_spline = csaps(f_param, f_tau, opt.smooth_tau);
            s_tau = fnval(tau_spline, s_param);
            val_spline = csaps(f_param, f_val, opt.smooth_val);
            s_val = fnval(val_spline, s_param);
            plot(hax, s_param, s_tau, 'k:', 'LineWidth', 1)
            
            coi = h2.ColorOrderIndex;
            loi = h2.LineStyleOrderIndex;
            plot(h2, f_param, f_tau, '.', 'MarkerSize', 8)
            h2.ColorOrderIndex = coi;
            h2.LineStyleOrderIndex = loi;
            plot(h2, s_param, s_tau, 'LineWidth', .5)
            
            plot(h3, f_param, f_val, '.', 'MarkerSize', 8)
            h3.ColorOrderIndex = coi;
            h3.LineStyleOrderIndex = loi;
            plot(h3, s_param, s_val, 'LineWidth', .5)
        end
    end
    
    % Plot x0 in black.
    if ~isempty(idx_x0)
        x0_feature = threads(idx_x0);
        param = x0_feature.param;
        x0 = x0_feature.val;
        s_param = (1:500)/500*(max(f_param) - min(f_param)) + min(f_param);     
        val_spline = csaps(param, x0, opt.smooth_val);
        s_val = fnval(val_spline, s_param);
        plot(param, x0, 'k.', 'MarkerSize', 8)
        plot(s_param, s_val, 'k-', 'LineWidth', .5)
    end
        
    txt1 = [self.name, '   Results of discretization'];
    if normalized
        txt1 = [txt1, '   ', '(normalized)'];
    end
    
    txt = {txt1, self.short_header};
    hax.Title.String = txt;
    hax.Title.FontSize = 9;
    hax.Title.FontWeight = 'normal';
    
    axes(h2)
    h2.Title.String = txt;
    h2.Title.FontSize = 9;
    h2.Title.FontWeight = 'normal';    
    xlb = opt.group_param_label;
    xlabel(xlb)
    ylabel(['\tau (', opt.time_unit, ')'])
    grid on
    
    axes(h3)
    line(xlim, [0, 0], 'Color', 'k', 'LineStyle', ':')
    xlabel(xlb)
    txt = 'Amplitude';
    if normalized
        txt = [txt, ' (nomalized)'];
    end
    ylabel(txt)
end

self.show_table('discrete')
end

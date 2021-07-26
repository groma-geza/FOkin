function diagnose(self)

history = self.history;
assert((history.DIAGNOSE > 0), ['No history is available. ',...
        'Run class function optimize() with options.DIAGNOSE '...
        '> 0 before calling diagnose().'])

if isfield(history, 'rho_numiter')
    iters = 1:self.numiter;
    hfig = self.h_diag_fig;
    if isempty(hfig)
        newfig = true;
    elseif ~ishandle(hfig)
        newfig = true;
    else
        newfig = false;
    end
    
    if newfig
        hfig = figure('Position', [100, 500, 1200, 500]);
        self.h_diag_fig = hfig;
        hax = axes(hfig);
    else
        figure(hfig)
        hax = findobj(hfig, 'Type', 'axes');
    end
    
    yyaxis left
    semilogy(hax, iters, history.r_norm, 'b-', iters, history.eps_pri,...
        'r-', iters, history.s_norm, 'b--', iters, history.eps_dual, 'r--')
    xlabel('Iteration')
    ylabel('r\_norm, eps\_pri, s\_norm, eps\_dual')
    grid on
    
    yyaxis right
    plot(hax, iters, history.nonzero, 'm')
    ylabel('Number of nonzero elements in z')
    
    hl = legend(hax, 'r\_norm', 'eps\_pri', 's\_norm', 'eps\_dual',...
        '# nonzero');
    hl.Location = 'northeast';
    
    hdc = datacursormode(hfig);
    hdc.UpdateFcn = @dcupdate;
    
    txt = [self.name, '   RHO0 = ', num2str(self.options.RHO0, '%.2E'),...
        '   iters = ', num2str(self.numiter)];
    hax.Title.String = txt;
    drawnow
    
    last = length(history.rho_numiter);
    if newfig
        first = 1;
        fprintf('\n')
        disp('History of rho-update:')
    else
        first = last;
    end
    for ii = first:last
        disp(['iter: ', num2str(history.rho_numiter(ii)),...
            '   ratio: ', num2str(history.rho_ratio(ii)),...
            '   rho: ', num2str(history.rho_value(ii))])
    end   
else
    fprintf('\n')
    disp('No new data for diagnosis')
end
end

%%
function txt = dcupdate(~, event_obj)
% Data cursor update function with high precision
pos = event_obj.Position;
txt = {sprintf('X: %.6g', pos(1)), sprintf('Y: %.6E', pos(2))};
end
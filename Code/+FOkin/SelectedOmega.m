classdef SelectedOmega < FOkin.SelectedModel
    %%
    properties (SetAccess = {?FOkin.SelectedModel, ?FOkin})
         best_omega,  best_value, best_omega_1STD
    end
    
    %%
    properties (SetAccess = {?FOkin.SelectedModel})
        omega_grid, average_support_size, step_right_limit
    end
        
    %%
    methods
        function count_support(self, ngrid, include_zero, verbose)
            arguments
                self
                ngrid(1,1) double {mustBeInteger,...
                    mustBeGreaterThan(ngrid, 2)} = 500;
                include_zero(1,1) logical = true;
                verbose(1,1) logical = false;
            end
                        
            fokin = self.fokin;
            options = fokin.options;
            fokin.optimizer.lambda = self.lambda_range;
            omega_range = self.omega_range;
                     
            lo = log10(omega_range(1));
            hi = log10(omega_range(2));
            if include_zero
                len = ngrid - 2;
            else
                len = ngrid - 1;
            end
            step = (hi - lo)/(len);
            log_omega_grid = lo:step:hi;
                       
            % Include self.best_omega in the grid
            log_best_omega = log10(self.best_omega);
            [~, idx] = min(abs(log_omega_grid - log_best_omega));
            log_omega_grid(idx) = log_best_omega;
            
            % Correct if any of the limiting values is lost.            
            if ~any(log_omega_grid == lo)
                log_omega_grid = [lo, log_omega_grid];
                ngrid = ngrid + 1;
                idx = idx + 1;
            end
            if ~any(log_omega_grid == hi)
                log_omega_grid = [log_omega_grid, hi];
                ngrid = ngrid + 1;
            end            
            omega_grid = 10.^log_omega_grid;                                
            
            if include_zero
                % Extent the grid with zero
                omega_grid = [0, omega_grid];
                idx = idx + 1;
            end
            
            average_support_size = zeros(1, ngrid);
            
            if verbose
                disp(['Running support size counting for ''',...
                    fokin.name, ''' in ' num2str(ngrid), ' steps...'])
                disp('Step in progress: ')
            end
            
            for jj = 1:ngrid              
                if verbose
                    fprintf(1, '%3.d ', jj)
                    if ((jj > 1) && (jj < ngrid) && ~mod(jj, 15))
                        fprintf(1, '\n')
                    end
                end
                fokin.optimizer.omega = omega_grid(jj); 
                fokin.optimize()
                res = fokin.results();
                support = (abs(res.x) > options.nonzero_limit);
                average_support_size(jj) = sum(support, 'all')/fokin.p;
            end
            
            if verbose
                fprintf(1, '\n')
            end
            self.omega_grid = omega_grid;
            self.average_support_size = average_support_size;
            
            level = self.average_support_size(idx);
            idx = find((average_support_size == level), 1, 'last');
            self.step_right_limit = omega_grid(idx);  
        end
        
        %%
        function show(self)
            show@FOkin.SelectedModel(self)
            if ~isempty(self.average_support_size)
                drawnow
                hfig = flip(findobj('Type', 'Figure')); % ordered by the 
                % time of creation
                figure(hfig(end-1)) % meanMSE vs. \omega
                hax = gca;
                yyaxis right
                color = [0.49, 0.18, 0.56]; % purple
                hax.YColor = color;
                % The first point of self.omega_grid is 0;
                hlin = semilogx(self.omega_grid(2:end),...
                    self.average_support_size(2:end), 'Color', color,...
                    'LineWidth', 2);
                hlin.DisplayName = 'Average support size';
                lims = ylim;
%                 lims(1) = 0;
                lims = [lims(1)-0.2, lims(end) + 0.2];
                ylim(lims)
                
                ylabel('Average support size')
                
                right = self.step_right_limit;
                hax.Title.String{2} = [hax.Title.String{2},...
                    '   \omega_{step right limit} = ', num2str(right,...
                    '%.2E')];
                hlin = line([right, right], ylim, 'Color', color,...
                    'LineStyle', ':',  'LineWidth', 1.5);
                hlin.DisplayName = 'Right limit of step selected';
            end
        end
        
    end    
end
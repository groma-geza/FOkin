classdef SelectedLambda < FOkin.SelectedModel
    %%
    properties (SetAccess = {?FOkin.SelectedModel, ?FOkin})
        best_lambda,  best_value, best_lambda_1STD
    end
    
    %%
    properties (SetAccess = {?FOkin.SelectedModel})
        lambda_grid, num_features, step_left_limit
    end
    
    %%
    methods
        function count_features(self, ngrid, verbose)
            arguments
                self
                ngrid(1,1) double {mustBeInteger,...
                    mustBeGreaterThan(ngrid, 1)} = 500;
                verbose(1,1) logical = false;
            end
                        
            fokin = self.fokin;
            fokin.optimizer.omega = self.omega_range;
            lambda_range = self.lambda_range;
                     
            lo = log10(lambda_range(1));
            hi = log10(lambda_range(2));
            step = (hi - lo)/(ngrid - 1);
            log_lambda_grid = lo:step:hi;

            % Include self.best_lambda in the grid
            log_best_lambda = log10(self.best_lambda);
            [~, idx] = min(abs(log_lambda_grid - log_best_lambda));
            log_lambda_grid(idx) = log_best_lambda;
            
            % Correct if any of the limiting values is lost.            
            if ~any(log_lambda_grid == lo)
                log_lambda_grid = [lo, log_lambda_grid];
                ngrid = ngrid + 1;
                idx = idx + 1;
            end
            if ~any(log_lambda_grid == hi)
                log_lambda_grid = [log_lambda_grid, hi];
                ngrid = ngrid + 1;
            end            
            lambda_grid = 10.^log_lambda_grid;              
            
            num_features = zeros(1, ngrid);
            
            if verbose
                disp(['Running feature counting for ''', fokin.name,...
                    ''' in ' num2str(ngrid), ' steps...'])
                disp('Step in progress: ')
            end
            
            for jj = 1:ngrid
                if verbose
                    fprintf(1, '%3.d ', jj)
                    if ((jj > 1) && (jj < ngrid) && ~mod(jj, 15))
                        fprintf(1, '\n')
                    end
                end
                fokin.optimizer.lambda = lambda_grid(jj); 
                fokin.optimize()
                res = fokin.results();
                disc = res.discretize();
                num_features(jj) = numel(disc.discrete.threads);
            end
            
            if verbose
                fprintf(1, '\n')
            end
            self.lambda_grid = lambda_grid;
            self.num_features = num_features;
            
            level = self.num_features(idx);
            idx = find((num_features == level), 1);
            self.step_left_limit = lambda_grid(idx);  
        end
        
        %%
        function show(self)
            show@FOkin.SelectedModel(self)
            if ~isempty(self.num_features)
                drawnow
                hfig = flip(findobj('Type', 'Figure')); % ordered by the 
                % time of creation
                figure(hfig(end-1)) % meanMSE vs. \lambda
                hax = gca;
                yyaxis right
                color = [0.49, 0.18, 0.56]; % purple
                hax.YColor = color;
                hlin = semilogx(self.lambda_grid, self.num_features, ...
                    'Color', color, 'LineWidth', 2);
                hlin.DisplayName = 'Number of features';
                ylabel('Number of features')
                lims1 = ylim;
                lims2 = [lims1(1)-0.2, lims1(end) + 0.2];
                ylim(lims2)
                hax.YTick = lims1(1):1:lims1(end);
                
                left = self.step_left_limit;
                hax.Title.String{2} = [hax.Title.String{2},...
                    '   \lambda_{step left limit} = ', num2str(left,...
                    '%.2E')];
                hlin = line([left, left], lims2, 'Color', color,...
                    'LineStyle', ':',  'LineWidth', 1.5);
                hlin.DisplayName = 'Left limit of step selected';
            end
        end
        
    end
end
classdef FOkinResults < matlab.mixin.Copyable
    % Results from the calculations by class FOkin
    %%
    properties (Access = {?FOkin.FOkinResults, ?FOkin})
        norm_orig
        optimizer_infotext1, optimizer_infotext2, optimizer_class
        level
    end
    
    %%
    properties
        % FOkin copies its own corresponding properties into these ones,
        % but they are publicly allowed to change. The original values
        % are also copied to the info property.
        name = '';
        options           
    end
    
    %%
    properties (SetAccess = {?FOkin.FOkinResults, ?FOkin})
        % results of optimization 
        operation, timestamp,time, group_param, tau, norm, fitdata,...
        lambda, omega, x, x0, numiter, runtime, runtime_total,...
        group_objval, MSE_total, info
    end
    
    %%
    methods (Access = {?FOkin.FOkinResults, ?FOkin})
        function self = FOkinResults(results)
            if nargin > 0
                % Make a copy of a superclass object
                assert(isa(self, class(results)),...
                    ['Input argument must be an object of ''',...
                    class(self), ''' or its any superclass.'])
                for prop = metaclass(results).PropertyList'
                    name = prop.Name;
                    if ~prop.Dependent
                        self.(name) = results.(name);
                    end
                end
            end
        end
    end
    
    %%
    methods
        function denorm(self)
            assert((self.level > 0), 'No results to denorm.')
            
            if (self.norm == 1)
                return
            end
            self.x = self.x*self.norm;
            self.x0 = self.x0*self.norm;
            self.MSE_total = self.MSE_total*(self.norm)^2;
            p = length(self.fitdata);
            for ii = 1:p
                self.fitdata(ii).b = self.fitdata(ii).b*self.norm;
                self.fitdata(ii).bfit = self.fitdata(ii).bfit*self.norm;
                self.fitdata(ii).fit = self.fitdata(ii).fit*self.norm;
                self.fitdata(ii).residual =...
                    self.fitdata(ii).residual*self.norm;
                self.fitdata(ii).MSE = self.fitdata(ii).MSE*(self.norm)^2;               
            end
            self.norm = 1;            
        end
        
        %%
        function renorm(self)
            assert((self.level > 0), 'No results to renorm.')
            if (self.norm == self.norm_orig)
                return
            end
            self.norm = self.norm_orig;               
            self.x = self.x/self.norm;
            self.x0 = self.x0/self.norm;
            self.MSE_total = self.MSE_total/(self.norm)^2;
            p = length(self.fitdata);
            for ii = 1:p
                self.fitdata(ii).b = self.fitdata(ii).b/self.norm;
                self.fitdata(ii).bfit = self.fitdata(ii).bfit/self.norm;
                self.fitdata(ii).fit = self.fitdata(ii).fit/self.norm;
                self.fitdata(ii).residual =...
                    self.fitdata(ii).residual/self.norm;
                self.fitdata(ii).MSE = self.fitdata(ii).MSE/(self.norm)^2;               
            end         
        end
        
        
        %%
        % Externally implemented
        show(self, selection)
        res = discretize(self)
        
    end

    %%    
    methods (Access = protected)

        function txt = short_header(self)

            if ~isempty(self.group_objval)
                objval = ['   objval = ',...
                    num2str(self.group_objval, '%.6E')];
                rntm1 = '   runtime = ';
            else
                objval = '';
                rntm1 = '   runtime_{total} = ';
            end
            
            runtime_ = self.runtime;
            if isempty(runtime_)
                rntm2 = '';
            else
                rntm2 = [num2str(sum(self.runtime), 5), ' s'];
            end

            txt = [self.optimizer_infotext1,...
                rntm1, rntm2, objval,...
                '   MSE_{total} = ', num2str(self.MSE_total, '%.4E')];
        end

        %%
        function txt = long_header(self, idx)

            header1 = [self.name, '   ',...
                self.optimizer_class, '   ',...
                self.options.group_param_name, ' = ',...
                num2str(self.group_param(idx)), ' ',...
                self.options.group_param_unit, '   ',...
                self.optimizer_infotext1];

            header2 = self.optimizer_infotext2;

            if (self.level == 2)
                if self.options.group_lasso
                    objval = self.group_objval;
                    numiter_ = self.numiter;
                    runtime_ = self.runtime;
                else
                    objval = self.fitdata(idx).objval;
                    numiter_ = self.numiter(idx);
                    runtime_ = self.runtime(idx);
                end
                unit = ' s';
            else
                objval = [];
                numiter_ = [];
                runtime_ = [];
                unit  = '';
            end
                

            header3 = ['numiter = ', num2str(numiter_),...
                '   runtime = ', num2str(runtime_, 5), unit,...
                '   objval = ', num2str(objval, '%.6E'),...
                '   MSE = ', num2str(self.fitdata(idx).MSE, '%.4E'),...
                '   MSE_{total} = ', num2str(self.MSE_total, '%.4E')];

            txt = {header1, header2, header3};
        end

        %%
        function hax = plottau2D(self, x)

            set(gcf, 'position', [25, 230, 620, 510])

            group_param_ = self.group_param;
            tau_ = self.tau;

            hax = axes;
            h = surface(group_param_, tau_, x);
            h.LineStyle = 'none';
            % h.FaceColor = 'interp';
            h.FaceColor = 'flat';
            colormap(jet(256))
            hax.YScale = 'log';
            hax.YLim = [tau_(1)/1.1, tau_(end)*1.1];
            maxval = max(abs(x), [], 'all');
            hax.CLim = [-0.15*maxval, 0.15*maxval];
            axis manual
            hax.TickDir = 'out';
            hax.Box = 'off';
            hax.SortMethod = 'childorder';
            xlabel(self.options.group_param_label);
            ylabel(['\tau (', self.options.time_unit, ')'])

            colorbar;
        end
        
        %%
        function self_copy = copyElement(self)            
        % Ensures deep copying.            
            self_copy = copyElement@matlab.mixin.Copyable(self);
            opt= self.options.copy();
            opt.owner = [];
            self_copy.options = opt;
            self_copy.fitdata= self.fitdata.copy();
        end
        
    end

end
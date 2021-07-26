classdef FOkinDiscretized < FOkin.FOkinResults
    % Result of FOkin.FOkinResults.discretize()
    properties (SetAccess = {?FOkin.FOkinResults})       
        discrete
    end

    %%
    methods (Access = {?FOkin.FOkinResults})     
        function self = FOkinDiscretized(results)
            arguments
                results(1,1) FOkin.FOkinResults
            end
            self@FOkin.FOkinResults(results)
        end
        
    end    
    %%
    methods
        %%
        function denorm(self)
            if (self.norm == 1)
                return
            end
            disc = self.discrete;
            limit = disc.neglected.limit;
            nmr = self.norm;
            if ~isempty(limit)
                self.reset() % temporarily for simplicity
            end
            
            for ii = 1:length(disc.threads)
                disc.threads(ii).val =...
                    disc.threads(ii).val*nmr;
            end
            
            if ~isempty(limit)
                self.neglect(limit)
            end
            disc.abs_max_amplitude =...
                disc.abs_max_amplitude*nmr;
            
            expfit = disc.expfit;
            if ~isempty(expfit)
                expfit.DADS = expfit.DADS*nmr;
                expfit.abs_max_amplitude =...
                    expfit.abs_max_amplitude*nmr;
                expfit.MSE_total = expfit.MSE_total*(nmr)^2;
                for ii = 1:length(expfit.fit)
                    expfit.fit{ii} = expfit.fit{ii}*nmr;
                    expfit.residual{ii} = expfit.residual{ii}*nmr;
                    expfit.MSE(ii) = expfit.MSE(ii)*(nmr)^2;
                end
            end
            
            denorm@FOkin.FOkinResults(self)
        end
        
        %%
        function renorm(self)
            if (self.norm == self.norm_orig)
                return
            end
            disc = self.discrete;
            limit = disc.neglected.limit;
            nmr_orig = self.norm_orig;
            if ~isempty(limit)
                self.reset() % temporarily for simplicity
            end
            
            for ii = 1:length(disc.threads)
                disc.threads(ii).val =...
                    disc.threads(ii).val/nmr_orig;
            end
            
            if ~isempty(limit)
                self.neglect(limit)
            end
            disc.abs_max_amplitude =...
                disc.abs_max_amplitude/nmr_orig;
            
            expfit = disc.expfit;
            if ~isempty(expfit)
                expfit.DADS = expfit.DADS/nmr_orig;
                expfit.abs_max_amplitude =...
                    expfit.abs_max_amplitude/nmr_orig;
                expfit.MSE_total = expfit.MSE_total/(nmr_orig)^2;
                for ii = 1:length(expfit.fit)
                    expfit.fit{ii} = expfit.fit{ii}/nmr_orig;
                    expfit.residual{ii} = expfit.residual{ii}/nmr_orig;
                    expfit.MSE(ii) = expfit.MSE(ii)/(nmr_orig)^2;
                end
            end
            
            renorm@FOkin.FOkinResults(self)
        end
        
        %%
        function neglect(self, limit)
            disc = self.discrete;
            threads = disc.threads;
            if isempty(threads)
                return
            end
            
            average_tau = disc.average_tau;
            rel_amplitude = disc.rel_amplitude;
            
            neglected = disc.neglected;
            
            % The elements of rel_amplitude are sorted in descending order
            k = length(threads);
            idx = [];
            for ii=1:k
                %     if (max(abs(threads{ii}.val)/maxval) <= limit)
                if rel_amplitude(ii) <= limit
                    idx = ii;
                    break
                end
            end
            
            if idx
                neglected.limit = limit;
                neglected.threads = [threads(idx:end), neglected.threads];
                neglected.average_tau = [average_tau(idx:end); neglected.average_tau];
                neglected.rel_amplitude = [rel_amplitude(idx:end);...
                    neglected.rel_amplitude];
                
                disc.threads = threads(1:(idx-1));
                disc.average_tau = average_tau(1:(idx-1));
                disc.rel_amplitude = rel_amplitude(1:(idx-1));
            end
            
            disc.neglected = neglected;
            
            self.service_expfit('neglect')            
        end
        
        %%
        function reset(self)            
            disc = self.discrete;
            threads = disc.threads;
            disc.threads = [threads, disc.neglected.threads];
            disc.neglected.threads = {};
            
            average_tau = disc.average_tau;
            disc.average_tau = [average_tau; disc.neglected.average_tau];
            disc.neglected.average_tau = [];
            
            rel_amplitude = disc.rel_amplitude;
            disc.rel_amplitude = [rel_amplitude; disc.neglected.rel_amplitude];
            disc.neglected.rel_amplitude = [];
            disc.neglected.limit = [];
            
            self.service_expfit('reset')
        end
        
        %%
        % Externally implemented
        do_expfit(self)
        showd(self)
        showexp(self, selection)
        
    end
    
    %%
    methods (Access = protected)
        %%
        function service_expfit(self, reason)
            threads = self.discrete.threads;
            expfit = self.discrete.expfit;
            if ~isempty(expfit.tau)
                if (length(threads) ~= length(expfit.tau))
                    self.discrete.expfit = FOkin.Expfit;
                    disp(['Results of exponential fitting is deleted ',...
                        'due to ', reason, '.'])
                end
            end
        end
        
        %%
        % Externally implemented
        show_table(self, type)
    end
   
end
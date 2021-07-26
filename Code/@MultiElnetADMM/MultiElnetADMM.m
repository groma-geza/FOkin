%%
classdef MultiElnetADMM < FOkin.MultiElnet
    %%
    properties (Access = protected)
        x, z, u, Atb
        numiter_, runtime_
        h_diag_fig
        rho
    end
    
    %%
    methods
        function self = MultiElnetADMM(varargin)
            self@FOkin.MultiElnet(varargin{:})
            self.options_ = MultiElnetADMMOptions;
            self.supports_group = true;
        end
        
        %%
        function numiter = numiter(self)
            numiter = self.numiter_;
        end
        
        %%
        function runtime = runtime(self)
            runtime = self.runtime_;
        end
        
        %%
        % Externally implemented
        diagnose(self)
    end
    
    %%
    methods (Access = {?FOkin.MultiElnet, ?FOkin, ?FOkin.FOkinResults})    
        %%
        function text = infotext2(self)
            opt = self.options;
            text = ['ABSTOL = ', num2str(opt.ABSTOL, '%.1E'),...
                '   RELTOL = ', num2str(opt.RELTOL, '%.1E'),...
                '   \rho0 = ', num2str(opt.RHO0, '%.1E'),...
                '   \alpha\_relax = ', num2str(opt.ALPHA_RELAX, '%.1f')];
        end
    end
    
    %%
    methods (Access = protected)
        
        function success = rebuild(self, varargin)
            success = rebuild@FOkin.MultiElnet(self, varargin{:});
            if success
                if self.multi_A
                    self.Atb = cell(1, self.p);
                    for ii = 1:self.p
                        self.Atb{ii} = (self.A{ii})'*self.b{ii};
                    end
                else
                    self.Atb = (self.A)'*self.b;
                end
            end
        end
 
        %%
        function result = do_get_result(self)
            result = self.z;
        end
        
        %%
        % Externally implemented
        do_reset(self)
        do_optimize(self)
        
    end
end


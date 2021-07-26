%%
classdef MultiElnetPDCO < FOkin.MultiElnet
       
    properties (Access = protected)
        x, y, z
        numiter_, runtime_
    end
    
    %%
    methods
        function self = MultiElnetPDCO(varargin)
            self@FOkin.MultiElnet(varargin{:})
            self.options = MultiElnetPDCOOptions;
            self.supports_group = false;
        end
        
        %%
        function numiter = numiter(self)
            numiter = self.numiter_;
        end
      
        %%
        function runtime = runtime(self)
            runtime = self.runtime_;
        end        
        
    end
    
    %%
    methods (Access = protected)
        function result = do_get_result(self)
            n = self.n;
            p = self.p;
            lngt = n*p;
            result = self.x(1:lngt) - self.x((lngt+1):(2*lngt));
            result = reshape(result, [n, p]);
        end
        
        %%
        % Externally implemented
        do_reset(self)
        do_optimize(self)
        
    end        
end


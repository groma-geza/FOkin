classdef SelectedModel < matlab.mixin.Copyable
    %%
    properties
        name
    end
    %%
    properties (SetAccess = {?FOkin.SelectedModel, ?FOkin})
        timestamp, lambda_range, omega_range, CV_spec, maxiter,...
            BO, runtime, info
    end
    
    %%
    properties (Access = {?FOkin.SelectedModel, ?FOkin})
        fokin
    end
    
    %%
    methods (Access = {?FOkin.SelectedModel, ?FOkin})
        function self = SelectedModel(selected)
            if nargin > 0
                % Make a copy of a superclass object
                assert(isa(self, class(selected)),...
                    ['Input argument must be an object of ''',...
                    class(self), ''' or its any superclass.'])
                for prop = metaclass(selected).PropertyList'
                    name = prop.Name;
                    if ~prop.Dependent
                        self.(name) = selected.(name);
                    end
                end
            end
        end
    end
    
    %%
    methods
        function show(self)
            if (~isscalar(self.lambda_range) &&...
                    ~isscalar(self.omega_range))
                show = 2;
            else
                show = 1;
            end 
            FOkin.process_selected(self, Inf, Inf, show);            
        end
    end
    
end
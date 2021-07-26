classdef Sections < matlab.mixin.Copyable
    %%
    properties (SetAccess = {?FOkin.SelectedLambdaOmega})
        omega_slice_width, best_lambda,  best_lambda_1STD
        lambda_slice_width, best_omega,  best_omega_1STD
        info
    end
    
    %%
    properties (Access = {?FOkin.SelectedLambdaOmega})
        parent
    end
    
    %%
    methods (Access = {?FOkin, ?FOkin.SelectedLambdaOmega})
        function self = Sections(parent, info)
            arguments
                parent FOkin.SelectedLambdaOmega
                info
            end
            self.parent = parent;
            self.info.previous = info;
        end
        
    end
        %%
    methods
        function show(self)
            FOkin.process_selected(self.parent, self.lambda_slice_width,...
                self.omega_slice_width, 1);
        end
        
    end

end
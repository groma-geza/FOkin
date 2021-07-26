classdef SelectedLambdaOmega < FOkin.SelectedModel
    %%
    properties (SetAccess = {?FOkin.SelectedModel, ?FOkin})
        best_lambda, best_omega, best_value, sections
    end
    
    %%
    methods
        function do_sections(self, lambda_slice_width,...
                omega_slice_width, show)
            arguments
                self
                lambda_slice_width(1,1) double {mustBePositive} = 0.4;
                omega_slice_width(1,1) double {mustBePositive} = 0.4;
                show(1,1) logical = false;
            end
            
            sections_ = FOkin.Sections(self, self.info);
            sections_.lambda_slice_width = lambda_slice_width;
            sections_.omega_slice_width = omega_slice_width;
            
            [sections_.best_lambda, sections_.best_omega,...
                sections_.best_lambda_1STD, sections_.best_omega_1STD] =...
                FOkin.process_selected(self, lambda_slice_width,...
                omega_slice_width, 0);
            self.sections = sections_;
            
            if show
                sections_.show()
            end
        end
        
    end
        
end
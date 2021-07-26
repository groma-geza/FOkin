classdef Fitdata < matlab.mixin.Copyable
    %%
    properties (SetAccess = {?FOkin.FOkinResults, ?FOkin})
        group_param, b, t0, fwhm, start, bfit, weightfit, fit, residual,...
            objval, MSE  
    end
end
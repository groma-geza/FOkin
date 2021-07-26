classdef Expfit < matlab.mixin.Copyable
    %%
    properties (SetAccess = {?FOkin.FOkinResults})
        tau, DADS, rel_amplitude, abs_max_amplitude,
        fit, residual, MSE, MSE_total
    end
end
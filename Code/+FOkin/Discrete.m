classdef Discrete < matlab.mixin.Copyable
    %%
    properties (SetAccess = {?FOkin.FOkinResults})
        threads, average_tau, rel_amplitude, abs_max_amplitude, neglected
        expfit
    end
end
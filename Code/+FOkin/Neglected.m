classdef Neglected < matlab.mixin.Copyable
    %%
    properties (SetAccess = {?FOkin.FOkinResults})
        threads, average_tau, rel_amplitude
        limit
    end
end
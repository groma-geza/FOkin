classdef MultiElnetADMMOptions < FOkin.MultiElnetOptions
    %%
    properties
        MAX_ITER(1,1) double {mustBePositive, mustBeInteger} = 1E5;
        ABSTOL(1,1) double {mustBePositive}   = 1E-10;
        RELTOL(1,1) double {mustBePositive}   = 1E-6;
        RHO0(1,1) double {mustBePositive} = 1.E-3;
        ALPHA_RELAX(1,1) double {mustBePositive} = 1.8;
    end
end
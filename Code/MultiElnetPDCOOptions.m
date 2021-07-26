classdef MultiElnetPDCOOptions < FOkin.MultiElnetOptions
    %%
    properties
        
        % values retutned by pdcoSet
%         MaxIter = 30;
%         FeaTol = 1.E-6;
%         OptTol = 1.E-6;
%         Print = 1;
%         StepTol = 0.99;
%         StepSame = 1;
%         x0min = 1;
%         z0min = 1;
%         mu0 = 1.E-1;
%         backtrack = 0;
%         Method = 0;
%         LSMRMaxIter = 10;
%         LSMRatol1 = 1e-10;
%         LSMRatol2 = 1e-15;
%         LSMRconlim = 1.0000e+12;
%         wait = 0;

        % Recommended values
        MaxIter(1,1) double {mustBePositive, mustBeInteger} = 1E4;
        FeaTol(1,1) double {mustBePositive} = 1.E-11;
        OptTol(1,1) double {mustBePositive} = 1.E-11;
        Print(1,1) logical = 0;
        StepTol(1,1) double {mustBeNonnegative,...
            mustBeLessThanOrEqual(StepTol,1)} = 1 - 1.E-4;
        StepSame(1,1) logical = 1;
        x0min(1,1) double {mustBePositive} = 0.001;
        z0min(1,1) double {mustBePositive} = 0.001;
        mu0(1,1) double {mustBePositive} = 1.E-12;
        backtrack(1,1) logical = 0;
        Method(1,1) double {mustBeMember(Method, 1)} = 1; % Sparse Cholesky
        LSMRMaxIter(1,1) double {mustBePositive} = 10;
        LSMRatol1(1,1) double {mustBePositive} = 1e-3;
        LSMRatol2(1,1) double {mustBePositive} = 1e-4;
        LSMRconlim(1,1) double {mustBePositive} = 1.0000e+12;
        wait(1,1) logical = 0;
               
    end
end
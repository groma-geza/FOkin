classdef RisefitOptions < matlab.mixin.Copyable    
    % Options for risefit
    %%
    properties
        individual(1,1) logical = false;
        separate_norm(1,1) logical = false;
        repeat(1,1) logical = false;
               
        elnet_lambda(1,1) double {mustBePositive} = 1.E-4;
        elnet_omega(1,1) double {mustBePositive} = 1.E-10;
        
        maxiter(1,1) double {mustBeInteger}= 100;
        
        show_details(1,1) {mustBeMember(show_details, [0, 1, 2, 3])} = 0;
        title_text char = '';
    end
end
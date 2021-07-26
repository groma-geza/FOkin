classdef MultiElnetOptions < matlab.mixin.Copyable
    %%
    properties
        DIAGNOSE(1,1) {mustBeMember(DIAGNOSE, [0, 1, 2])} = 0;
    end
end
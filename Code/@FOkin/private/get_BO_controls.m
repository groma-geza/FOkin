function controls = get_BO_controls(options, maxiter)
% Get Name, Value input parameter pairs for bayesopt() from the properties
% of options with the name of a leading 'bo_' string.
% input:  option   - an object of class FOkinOptions
%         maxiter  - value of 'MaxObjectiveEvaluations' in controls
% output: controls - a vector of cells array, containing Name, Value pairs
arguments
    options(1,1) FOkinOptions
    maxiter(1,1) double {mustBePositive, mustBeInteger}
end
mc = metaclass(options);
props = mc.PropertyList;
n = numel(props); 
controls = {};
for ii = 1:n
    get_access = props(ii).GetAccess;
    if ischar(get_access)
        if strcmp(get_access, 'public')
            name = props(ii).Name;
            value = options.(name);
            include = false;
            if startsWith(name, 'bo_')
                if ~isnumeric(value)
                    include = true;
                elseif ~any(isnan(value))
                    include = true;
                end
            end
            if include
                controls = [controls, {name(4:end)}, {value}]; %#ok<AGROW>
            end
        end
    end
end
controls = [controls, {'MaxObjectiveEvaluations'}, {maxiter}];
end
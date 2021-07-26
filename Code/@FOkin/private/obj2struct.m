function struct_obj = obj2struct(obj)
% Make a copy of a class object obj into a nested structure struct_obj.
% Only properties of public GetAccess are copied.
% Properties with a value of struct type are not supported.
struct_obj = [];
mc = metaclass(obj);
for prop = mc.PropertyList'
    if strcmp(prop.GetAccess, 'public')
        name = prop.Name;
        val = obj.(name);
        if isa(val, 'struct')
            error(['Properties with a value of struc type are not ',...
                'supported.'])
        end
        mc2 = metaclass(val);
        if ~isempty(mc2.PropertyList)
            val = FOkin.obj2struct(val);
        end
        struct_obj.(prop.Name) = val;
    end
end
end




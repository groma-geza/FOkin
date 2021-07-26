classdef MultiElnet < matlab.mixin.Copyable
    %%
    properties
        name char = '';
        A
        b
        % control parameters
        lambda(1,1) double {mustBePositive} = 1.E-4;
        alpha(1,1) double {mustBeNonnegative,...
            mustBeLessThanOrEqual(alpha,1)} = 1;

    end

    %%
    properties (Dependent = true)
        omega
    end
    %%
    properties (SetAccess = protected, Dependent = true)
        lambda1, lambda2 %control parameters       
    end
    
    %%
    properties (Dependent = true)
        options
    end
    
   %% 
    properties (SetAccess = protected)        
        supports_group
    end

    %%
    properties (SetAccess = protected, Dependent = true)
        result % result of the elastic net optimization
    end
    
    %%
    properties (Access = protected)
        n, m, p, multi_A
        old_n, old_p
        options_keys
        optimized
        history  % structure with parameters of optimization history
                 % defined in the subclasses
                 % and controlled by self.options.DIAGNOSE
        options_
    end
        
    %%
    methods
        function self = MultiElnet(A, b)
            arguments
                A = [];
                b = [];
            end
            
            self.A = A;
            self.b = b;
            
            % default parameter values
            self.old_n = [];
            self.old_p = [];
            self.options_ = FOkin.MultiElnetOptions();           
            self.history = struct();
            self.optimized = false;
            self.supports_group = false;
        end
            
        %%        
        function omega = get.omega(self)
            omega = 1 - self.alpha;
        end
        
        %%
        function set.omega(self, omega)
            arguments
                self
                omega(1,1) double {mustBeNonnegative,...
                    mustBeLessThanOrEqual(omega,1)}
            end
            self.alpha = 1 - omega;
        end
        
        %%
        function options = get.options(self)
            options = self.options_;
        end
        
        %%
        function set.options(self, options)
            classname = class(self.options_);
            assert(isa(options, classname), ['The ''options'' propery',...
                ' must be a ', classname, ' object.'])
            self.options_ = options;
        end
        
        %%
        function lambda1 = get.lambda1(self)
            lambda1 = self.lambda*self.alpha;
        end
           
        %%
        function lambda2 = get.lambda2(self)
            lambda2 = self.lambda*(1 - self.alpha);
        end
            
        %%
        function result = get.result(self)
            % Silly trick to ensure overwriting get.result in subclasses.
            if self.optimized
                result = self.do_get_result();
            else
                result = [];
            end
        end
        
        %%
        function optimize(self, warm)
            arguments
                self
                warm(1,1) logical = false;
            end
            if ~self.rebuild(warm)
                error('Problem rebuilding was not successful')
            end
            self.do_optimize
            self.optimized = true;
        end
        
        %%
        function reset(self)
            self.optimized = false; % no result is available
            % Indicate a need to execute self.do_reset() in self.rebuild()
            self.old_n = [];
            self.old_p = [];
        end
    end
    methods (Access = {?FOkin.MultiElnet, ?FOkin, ?FOkin.FOkinResults})
        %%
        function text = infotext1(self)
            text = ['\lambda = ', num2str(self.lambda,...
                '%.2E'),...
                '   \omega = ', num2str(self.omega, '%.2E')];
        end
        
        %%
        function text = infotext2(~)
            text = '';
        end
        
    end
    
    %%
    methods (Access = {?FOkin.MultiElnet, ?FOkin})
        %%
        % Externally implemented
        objval = objective(self, x)
    end
    
    %%
    methods (Access = protected)
        
        function self_copy = copyElement(self)            
        % Ensures deep copying.
            self_copy = copyElement@matlab.mixin.Copyable(self);
            self_copy.options_ = self.options_.copy();
        end
        
        %%
        function success = rebuild(self, warm)
            arguments
                self
                warm(1,1) logical = false;
            end
            A = self.A;
            b = self.b;
            self.optimized = false;
            if isempty(A) || isempty(b)
                success = false;
                return
            end
            incompatible = '''A'' and ''b'' have incompatible size.';
            if iscell(A)
                self.multi_A = true;
                if ~iscell(b)
                    error('''A'' is a cell but ''b'' is not.')
                end   
                if any(cellfun(@isempty, A))
                    error('Empty object(s)in cell ''A''.')
                end
                if any(cellfun(@isempty, b))
                    error('Empty object(s)in cell ''b''.')
                end                                
                p = length(A);
                if length(b) ~= p
                    error(incompatible)
                end
                self.p = p;
                self.m = zeros(1, p);
                self.n = size(A{1}, 2);
                for ii = 1:p
                    [self.m(ii), n] = size(A{ii});
                    if (n ~= self.n)
                        error(['Matrices in cell ''A'' have different ',...
                            'number of columns.'])
                    end
                    if size(b{ii}, 1) ~= self.m(ii)
                        error(incompatible)
                    end
                end
            else
                % Suppose that A is a matrix and keep in that form.
                % In this case b also must be a matrix.
                self.multi_A = false;
                [self.m, self.n] = size(A);
                if iscell(b)
                    error('''b'' is a cell but ''A'' is not.')
                end      
                if isempty(A)
                    error('Matrix ''A'' is an empty object.')
                end
                if isempty(b)
                    error('Vector ''b'' is an empty object.')
                end      
                [m, self.p] = size(b);
                if (m ~= self.m)
                    error(incompatible)
                end
            end
                                 
            % execute cold restart if needed
            if (isempty(self.old_n) || isempty(self.old_p))
                % the method is executed for the first time
                % or after a self.reset()
                warm = false;
            end
            if ~warm
                self.old_n = self.n;
                self.old_p = self.p;
                self.do_reset()
            elseif (self.old_n ~= self.n) || (self.old_p ~= self.p)
                error(['Warm rebuilding is not possible due to ',...
                    'size change.'])
            end
            
            success = true;
        end
        
    end
    
%%
    methods (Abstract = true)
        numiter(self)
        runtime(self)
    end

%%
    methods (Abstract = true, Access = protected)
        do_reset(self) % sets internal values for cold restart
        result = do_get_result(self)
        do_optimize(self)
    end
    
end


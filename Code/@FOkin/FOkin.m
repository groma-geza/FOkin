
classdef FOkin < matlab.mixin.Copyable
    %%
    properties
    end
    
    %%
    properties (SetAccess = protected)
        data double
        time double
        group_param double
    end
    
    %%
    properties (Dependent = true)
        weight, t0, fwhm, start, name
        options, optimizer
    end
    
    %%
    properties (SetAccess = protected)      
        tau, A, b, bfit, 
        n, m, p
    end
    
    %%
    properties (Access = protected)
        data_norm, norm, multi_A
        ready, numiter, runtime
        weight_ 
        t0_ 
        fwhm_
        start_ 
        name_
        options_
        optimizer_
        x_norm, objval
        weightfit
        results_
        select_model_running
    end
    
    %%
    methods
        function  self = FOkin(varargin)
            self.recreate(varargin{:})
        end
        
        %%
        function weight = get.weight(self)
            weight = self.weight_;
        end
        
        %%
        function set.weight(self, weight)
            arguments
                self
                weight double {mustBePositive,...
                     mustBeLessThanOrEqual(weight,1)}
            end
            m = self.m;
            p = self.p;
            if isscalar(weight)  % mxp matrix
                weight = ones(m, p)*weight;
            end
            [row, col] = size(weight);

            assert((row == m), ['Argument ''weight'' and ',...
                '''self.data'' have incompatible size.'])
            if col == 1
                weight = repmat(weight, 1, p);
            elseif col ~= p
                error(['Argument ''weight'' and ''self.data'' have ',...
                    'incompatible size.'])
            end
                
            self.weight_ = weight;
            self.rebuild()
        end
        
        %%
        function t0 = get.t0(self)
            t0 = self.t0_;
        end
        
        %%
        function set.t0(self, t0)
            arguments
                self
                t0 double {mustBeNonnegative} 
            end
            if ~isscalar(t0)
                row = size(t0, 1);
                if row > 1
                    t0 = t0';
                end
                [row, col] = size(t0);
                assert((row == 1), 'Argument ''t0'' must be a vector.')
                assert((col == self.p), ['Argument ''t0'' and ',...
                    '''self.group_param'' have incompatible size.'])
            end
            self.t0_ = t0;
            self.rebuild()
        end
        
        %%
        function fwhm = get.fwhm(self)
            fwhm = self.fwhm_;
        end
        
        %%
        function set.fwhm(self, fwhm)
            arguments
                self
                fwhm double {mustBeNonnegative}
            end
            if ~isscalar(fwhm)
                row = size(fwhm, 1);
                if row > 1
                    fwhm = fwhm';
                end
                [row, col] = size(fwhm);

                assert((row == 1), 'Argument ''fwhm'' must be a vector.')
                assert((col == self.p), ['Arguments ''fwhm'' and ',...
                    '''self.group_param'' have incompatible size.'])
            end
            self.fwhm_ = fwhm;
            self.rebuild()
        end

        %%
        function start = get.start(self)
            start = self.start_;
        end
        
        %%
        function set.start(self, start)
            arguments
                self
                start double {mustBeInteger, mustBeNonnegative}
            end
            if ~isscalar(start)
                row = size(start, 1);
                if row > 1
                    start = start';
                end
                [row, col] = size(start);
                assert((row == 1), 'Argument ''start'' must be a vector.')
                assert((col == self.p), ['Arguments ''start'' and ',...
                    '''self.group_param'' have incompatible size.'])
            end
            self.start_ = start;
            self.rebuild()
        end
        
        %%
        function name = get.name(self)
            name = self.name_;
        end
        
        %%
        function set.name(self, name)
            arguments
                self
                name char
            end
            self.name_ = name;
            self.rebuild()
        end
        
        %%
        function options = get.options(self)
            options = self.options_;
        end
        
        %%
        function set.options(self, options)
            arguments
                self
                options(1,1) FOkinOptions
            end
            options.owner = self;
            self.options_ = options;
            self.rebuild()
        end
        
        %%
        function optimizer = get.optimizer(self)
            optimizer = self.optimizer_;
        end
        
        %%
        function set.optimizer(self, optimizer)
            arguments
                self
                optimizer(1,1) FOkin.MultiElnet
            end
            self.optimizer_ = optimizer;
            self.rebuild()
        end
        
        %%
        function res = results(self, info)
            arguments
                self
                info(1,1) struct = struct;
            end
            % Get results of optimization.
            
            if isempty(self.select_model_running)
                self.select_model_running = false;
            end
            
            if (isempty(self.results_) || self.select_model_running)
                self.build_results('no operation')
                self.select_model_running = false;
            end
            
            res = self.results_;
            res.info.previous = info;
        end
                         
        % Externally implemented
        optimize(self, varargin)
        
        [meanMSPE, constraints, user_data] = CV(self, lambda, omega,...
            method, warm)

        selected = select_model(self, lambda, omega, maxiter, info,...
            method, warm)
        
        [t0, fwhm, rel_params] = calc_t0_fwhm(self, maxiter)
    end
 
    %%
    methods (Access = protected)
        %%
        function  recreate(self, data, time, group_param, weight, t0,...
                fwhm, start, name)
            arguments
                self
                data double
                time double {mustBePositive}
                group_param double
                weight double {mustBePositive,...
                     mustBeLessThanOrEqual(weight,1)} = 1;
                t0 double {mustBeNonnegative} = 0;
                fwhm double {mustBeNonnegative} = 0;
                start double {mustBePositive, mustBeInteger}= 1;
                name char = '';
            end
            
            % Input
            %  data       - mxp matrix, containing the kinetic data in
            %               arbitrary units. A row of the matrix
            %               corresponds to a given group_param value and
            %               different time values and a column corresponds
            %               to a  given time value and different
            %               group_param values.
            % time        - mx1 vector, containing the time points in
            %               increasing order but in arbitrary
            %               steps. An arrangement close to logarithmically
            %               equidistant point is suggested.
            %               group_param - 1xp vector, containing a group
            %               parameter (typically wavelength) in increasing
            %               order but in arbitrary steps.
            %
            % weight      - scalar, mx1 vector or
            %               mxp matrix, containing the weight of fitting at
            %               different points of time. A vector means equal
            %               weights for all columns of data, a scalar means
            %               equal values for all data points. Value of 1
            %               means unweighted fitting. Default value is 1.
            % t0          - scalar or 1xp vector, containing the the real
            %               t=0 values (on the timescale of time),
            %               at the different points of group_param. A
            %               scalar means equal values for all points.
            %               Default value is 0.
            % fwhm        - scalar or 1xp vecto1, containing the value of
            %               the FWHM of the Gaussian describing
            %               the temporal IRF of the measuring device at the
            %               different points of group_param. A scalar means
            %               equal values for all points. Zero can be
            %               specified for an instantaneous IRF. Default
            %               value is 0.
            % start       - scalar or 1xp vector, containing the index of
            %               time from where the fitting starts at
            %               the different points of group_param.  A scalar
            %               means equal values for all points. Default
            %               value is 1.
            % name        - char array describing the dataset.
            
            self.data = data;
            self.time = time;
            self.group_param = group_param;
            
            % Check and adjust the dimensions of the arguments.
            [m, col] = size(self.time); % column vector
            if col > 1
                if m > 1
                    error(['Argument ''time has'' to be a vector not',...
                        ' a matrix.'])
                else
                    self.time = self.time';
                    m = col;
                end
            end
            self.m = m;
            
            [row, p] = size(self.group_param); % row vector
            if row > 1
                if p > 1
                    error(['Argument ''group_param'' has to be a ',...
                        'vector not a matrix.'])
                else
                    self.group_param = self.group_param';
                    p = row;
                end
            end
            self.p = p;
            
            [row, col] = size(self.data); % mxp matrix
            if (row ~= m)|| (col ~= p)
                self.data = self.data';
                [row, col] = size(self.data);
                assert((row == m), ['Arguments ''data'' and ''time'' '...
                    'have incompatible size.'])
        
                assert((col == p), ['Arguments ''data'' and',...
                    '''group_param'' have incompatible size.'])
            end
            
            % Store normalized data in self.data_norm
            self.norm = max(max(abs(self.data)));
            self.data_norm = self.data/self.norm;
            
            % weight, t0, fwhm and start are dependent on the previously
            % define properties          
            self.weight = weight;
            self.t0 = t0;
            self.fwhm = fwhm;
            self.start = start;
            self.name = name;
            
            opt = FOkinOptions();
            opt.owner = self;
            self.options_ = opt;
            self.optimizer_ = MultiElnetADMM();
            self.rebuild()
        end                

        %%
        function self_copy = copyElement(self)            
        % Ensures deep copying.            
            self_copy = copyElement@matlab.mixin.Copyable(self);
            opt = self.options_.copy();
            opt.owner = self_copy;
            opt.loadobj(opt);
            self_copy.options_ = opt;
            self_copy.optimizer_ = self.optimizer_.copy();
            self_copy.results_ = self.results_.copy();
        end
        
        %%
        % Externally implemented
        rebuild(self)
        build_results(self, operation)
    end
        
end



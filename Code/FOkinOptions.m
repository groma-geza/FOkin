classdef FOkinOptions < matlab.mixin.Copyable
    % Options for FOkin and related classes
    
    %%
    properties (Access = {?FOkin, ?FOkin.FOkinResults})
        owner = [];
    end
    %%
    properties (SetObservable, AbortSet)
        % Invokes the rebuild() method of the FOkin class via an event
        % listener if any of these properties is changed.
        n_tau_decade = 50; % Number of points in a decade of the tau
                           % property of the FOkin object.
        extension_lo = 1; % Number of decades the lower limit of tau
                          % will be extended with respect to the first point
                          % of the time property of the FOkin object.
                          % Fractional and negative values are allowed.
        extension_hi = 1; % Number of decades the higher limit of tau
                          % will be extended with respect to the last point
                          % of the time property of the FOkin object.
                          % Fractional and negative values are allowed.    
    end
    
    %%
    properties
        %%
        % General options
        num_par_workers = 0; % Number of required parpool workers in
                             % parfor loops. Set to zero to run parfor as
                             % for. Set to Inf to use all the available
                             % workers. Nonzero numbers need the Parallel
                             % Computing Toolbox of MATLAB to be installed.
        parpool_spec = 'local'; % Ignored if num_par_workers == 0.
                             % Otherwise, a char array representing a valid
                             % profile name acceptable as the 'resources'
                             % input argument of the parpool() function in
                             % the Parallel Computing Toolbox.
        group_lasso = true; % Applies only if the 'p' property of the 
                            % parent FOkin object is > 1. If true the whole
                            % dataset is fitted together, applying the
                            % group-lasso penalty in the objective
                            % function. If false the different subsets of
                            % data (colums of 'bfit' if it is a vector or
                            % the vectors stored in its elements if it is a
                            % cell array)are fitted separately with the
                            % lasso penalty.
        signal_label = 'Fluorescence (rel)'; % used in plot y labels
        time_unit = 'ps';   % used in plot x labels and table headers
        group_param_label = 'Wavelength (nm)'; % used in plot x labels
        group_param_name = 'wavelength'; % used in plot titles
        group_param_unit = 'nm'; % used in plot titles

        % Options for cross-validation
        cv_nfold = 10; % Number of folds in k-fold CV.
        cv_nrep = 1E4; % Number of repetitions in RCVnv.
        nc_rel = 0.9; % n_c/n where n_c is the sample size of the
                      % construction subset and n is the total sample size
                      % in RCVnv;
        lsqminnorm_tol = []; % Optional tolerance for the lsqminnorm()
                             % function of MATLAB applied in the RCVnv
                             % algorithm. If empty (recommended) the
                             % tolerance is determined internally.
        nonzero_limit = 0; % An element of the solution matrix X is 
                           % considered nonzero if its absolute value is
                           % higher than this value. Applied in the
                           % following methods: FOkin.CV (for RCVnv),
                           % FOkin.FOkinResults.discretize and
                           % FOkin.SelectedOmega.count_support. Suggested
                           % value: 0 for MultiElnetADMM
                           %        1.E-6 for MultiElnetPDCO.
        
        %%
        % Options for FOkin.FOkinResults.discretize()
        
        min_feature_gap = 1 % Minimum number of points in a gap required
                            % for separating two features.
        thread_row_region = 20; % Maximum number of rows to step up or
                                % down for finding the next point of a
                                % thread.
        min_thread_length = 3; % Minimum allowed length of a thread.
        
        thread_colors = [0   0    1      % blue
            0.39 0.83 0.07   % green
            1    0    0      % red
            0.06 1    1      % cyan
            1    0    1      % magenta
            0.93 0.69 0.13   % gold
            0.64 0.08 0.18]; % brown
            % color order used in FOkin.FOkinDiscretized.showd()
        
        % Smoothness of the splines defining the continuous line across the
        % time constants and the amplitudes, respectively, on a scale of
        % (0-1). 0: a constant equal to the mean of the values 
        %        1: the smoothest spline crossing all points.
        smooth_tau = 1.E-5;
        smooth_val = 1.E-3;
        
        %%
        % Options for controlling the execution of the bayesopt() function
        % in the Statistics and Machine Learning Toolbox of MATLAB, applied
        % in the FOkin.select_params() and FOkin.calc_t0_fwhm() methods.
        % The optional Name, Value input argument pairs of bayesopt() are
        % formed from the property names (without the leading 'bo_') and
        % their values, respectively. See the corresponding MATLAB
        % documentation for the details. A NaN value forces of using the
        % default MATLAB value. 'MaxObjectiveEvaluations' is not included
        % here as it obtains its value directly from an input argument of
        % the methods.
        
        bo_AcquisitionFunctionName = 'expected-improvement-plus';
        bo_IsObjectiveDeterministic = false;
        bo_ExplorationRatio = 0.5;
        bo_GPActiveSetSize = 1E3; % slower but more precise optimization
        bo_UseParallel = false;
        bo_ParallelMethod = 'clipped-model-prediction';
        bo_MinWorkerUtilization = NaN;
        bo_MaxTime = Inf;
        bo_NumSeedPoints = 4;
        bo_XConstraintFcn = [];
        bo_ConditionalVariableFcn = [];
        bo_NumCoupledConstraints = 0;
        bo_AreCoupledConstraintsDeterministic = NaN;
        bo_CoupledConstraintTolerances = NaN;
        bo_Verbose = 1;
        bo_OutputFcn = {};
        bo_SaveFileName = 'BayesoptResults.mat';
        bo_SaveVariableName = 'BayesoptResults';
        bo_PlotFcn = {@FOkin.FOkinPlotObjectiveModel,...
            @plotAcquisitionFunction};
        bo_InitialX = NaN;
        bo_InitialObjective = [];
        bo_InitialConstraintViolations = [];
        bo_InitialErrorValues = [];
        bo_InitialUserData = [];
        bo_InitialObjectiveEvaluationTimes = [];
        bo_InitialIterationTimes = {};
        
        %%
        % Options for FOkin.calc_t0_fwhm()
        n_t0_knots = 3;         % Number of knots in the spline
                                % defining the values of output argument
                                % vector t0.
        t0_range = [0.5, 1];    % Allowed range of t0 relatively to its 
                                % automatically estimated value.
        n_fwhm_knots = 1;       % Number of knots in the spline
                                % defining the values of output argument
                                % vector fwhm. Applies only if fwhm_fixed
                                % is empty.
        fwhm_range = [35, 50];  % Allowed range of fwhm in its own units. 
                                % Applies only if fwhm_fixed is empty.
        fwhm_fixed = [];        % Vector of the fixed values of fwhm.
                                % If empty fwhm will be approximated by
                                % a spline.
    end
    
    %%
    methods
        function self = FOkinOptions
            self.loadobj(self);  % activates event listener
        end
    end
    %%
    methods (Access = protected)
        function update_owner(self, ~, ~)
            % called only by event listener
            if ~isempty(self.owner)
                self.owner.options = self;
            end
        end
    end
    
    %%
    methods (Static)
        function self = loadobj(self)
            % Activates event listener on loading of a class object.
            % Called also in the constructor
            mc = metaclass(self);
            proplist = mc.PropertyList;
            proplist = proplist([proplist.SetObservable]);
            addlistener(self, proplist, 'PostSet', @self.update_owner);
        end
    end
    
end
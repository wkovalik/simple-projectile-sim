classdef SensorModel < handle
    properties
        sensorID
        params

        noiseCovar

        projectileModel
        earthModel

        estimatedParams
        estimatedParamSourceMap
        estimatedParamAnalyticJacobianMap

        jacobianMethod
    end

    properties (Dependent)
        invNoiseCovar

        N_PROJECTILE_STATES
        N_ESTIMATED_PARAMS
    end

    properties (Abstract, Constant)
        N_MEASUREMENTS
    end

    properties (Constant)
        PARAM_SOURCE_MAP = dictionary( ...
                                 "CD",     "projectileModel", ...
                                 "vWindx", "earthModel",      ...
                                 "vWindy", "earthModel"       ...
                             );
        PARAM_ANALYTIC_JACOBIAN_MAP = dictionary( ...
                                          "CD",     @computeAnalyticDragJacobian,  ...
                                          "vWindx", @computeAnalyticWindxJacobian, ...
                                          "vWindy", @computeAnalyticWindyJacobian  ...
                                      );
    end

    methods
        % Constructor ==============================================================================
        function self = SensorModel()
            self.noiseCovar = zeros(self.N_MEASUREMENTS);
            
            % Initialize list of parameters to be estimated (used to compute parameter Jacobian)
            self.estimatedParams = [];

            % Set initial Jacobian computation method
            self.jacobianMethod = Constants.DEFAULT_JACOBIAN_METHOD;
        end


        % Jacobian methods =========================================================================
        function H = computeNumericStateJacobian(self, state)
            H = zeros(self.N_MEASUREMENTS, self.N_PROJECTILE_STATES);
            
            % Numerically compute partial derivatives w.r.t. each projectile state
            for i = 1:self.N_PROJECTILE_STATES
                % Compute state perturbation
                delta = Constants.JACOBIAN_PERTURBATION_FACTOR * (1 + abs(state(i)));
                
                % Compute forward-perturbed measurement
                statePlus = state;
                statePlus(i) = statePlus(i) + delta;

                measurementPlus = self.computeMeasurement(statePlus);
                
                % Compute backward-perturbed measurement
                stateMinus = state;
                stateMinus(i) = stateMinus(i) - delta;

                measurementMinus = self.computeMeasurement(stateMinus);
                
                % Approximate partial derivative using central difference method
                H(:, i) = (measurementPlus - measurementMinus) / (2 * delta);
            end
        end

        function H = computeNumericParamJacobian(self, state)
            H = zeros(self.N_MEASUREMENTS, self.N_ESTIMATED_PARAMS);

            % Numerically compute partial derivatives w.r.t. each estimated parameter
            for i = 1:self.N_ESTIMATED_PARAMS
                paramName = self.estimatedParams{i};
                paramSource = self.estimatedParamSourceMap{i};

                param = self.(paramSource).params.(paramName).value;

                % Compute parameter perturbation
                delta = Constants.JACOBIAN_PERTURBATION_FACTOR * (1 + param);

                % Compute forward-perturbed measurement
                paramPlus = param + delta;
                self.(paramSource).params.(paramName).value = paramPlus;

                measurementPlus = self.computeMeasurement(state);

                % Compute backward-perturbed measurement
                paramMinus = param - delta;
                self.(paramSource).params.(paramName).value = paramMinus;

                measurementMinus = self.computeMeasurement(state);

                % Approximate partial derivative using central difference method
                H(:, i) = (measurementPlus - measurementMinus) / (2 * delta);

                % Reset parameter
                self.(paramSource).params.(paramName).value = param;
            end
        end

        function H = computeAnalyticParamJacobian(self, state)
            H = zeros(self.N_MEASUREMENTS, self.N_ESTIMATED_PARAMS);
            
            % Build Jacobian using individual Jacobians corresponding to each estimated parameter
            for i = 1:self.N_ESTIMATED_PARAMS
                H(:, i) = self.estimatedParamAnalyticJacobianMap{i}(self, state);  % See Note 1
            end
        end


        % Getters ==================================================================================
        function invNoiseCovar = get.invNoiseCovar(self)
            invNoiseCovar = inv(self.noiseCovar);
        end

        function N_PROJECTILE_STATES = get.N_PROJECTILE_STATES(self)
            N_PROJECTILE_STATES = self.projectileModel.state.N_STATES;
        end

        function N_ESTIMATED_PARAMS = get.N_ESTIMATED_PARAMS(self)
            N_ESTIMATED_PARAMS = length(self.estimatedParams);
        end


        % Setters ==================================================================================
        % TODO: set.params validation

        function set.projectileModel(self, projectileModel)
            self.projectileModel = Validator.validateType(projectileModel, "Projectile");
        end

        function set.earthModel(self, earthModel)
            self.earthModel = Validator.validateType(earthModel, "Earth");
        end

        function set.estimatedParams(self, estimatedParams)
            if ~isempty(estimatedParams)
                self.estimatedParams = Validator.validateType(estimatedParams, "string");
            end
            
            % Rebuild mapping lists for all estimated parameters
            self.buildEstimatedParamMaps();
        end

        function buildEstimatedParamMaps(self)
            % Builds list of object source location and list of Jacobian compute function handles
            % for each corresponding estimated parameter. See Note 1
            % Example: self.estimatedParams = [ "CD", "vWindx" ]
            %          self.estimatedParamSourceMap = { "projectile", "earth" }
            %          self.estimatedParamAnalyticJacobianMap = { @self.computeAnalyticDragJacobian, @self.computeAnalyticWindxJacobian }

            self.estimatedParamSourceMap = {};
            self.estimatedParamAnalyticJacobianMap = {};

            if ~isempty(self.estimatedParams)
                for i = 1:self.N_ESTIMATED_PARAMS
                    estimatedParam = self.estimatedParams(i);
                    
                    self.estimatedParamSourceMap{end + 1} = self.PARAM_SOURCE_MAP(estimatedParam);
                    self.estimatedParamAnalyticJacobianMap{end + 1} = self.PARAM_ANALYTIC_JACOBIAN_MAP(estimatedParam);
                end
            end
        end
    end

    methods (Abstract)
        y = computeMeasurement(self, state)      
        H = computeAnalyticStateJacobian(self, state)
    end
end


% Note 1
%
% The original implementation of the computeAnalyticParamJacobian() function was
%
% | for i = 1:self.N_ESTIMATED_PARAMS
% |     estimatedParam = self.estimatedParams(i);
% | 
% |     A(:, i) = feval(self.PARAM_ANALYTIC_JACOBIAN_MAP(estimatedParam), state);
% | end
%
% However, directly looking up the Jacobian function handles in self.PARAM_ANALYTIC_JACOBIAN_MAP
% ended up being quite slow (around twice as slow as the current implementation, in fact). I dunno
% why. My best guess is the local, individual self.computeAnalyticParamJacobian function handles
% were being rebuilt on every single lookup of self.PARAM_ANALYTIC_JACOBIAN_MAP(estimatedParam).
%
% So instead, a local self.estimatedParamAnalyticJacobianMap property is built, which is a cell
% array of the subset of Jacobian function handles in the full self.PARAM_ANALYTIC_JACOBIAN_MAP for
% the parameters in self.estimatedParams only. This cell array lists the handles in the correct
% order according to the order in self.estimatedParams and is rebuilt whenever self.estimatedParams
% is updated.
%
% For example:
% 
% | self.estimatedParams = ["CD", "vWindx"]
% | self.estimatedParamAnalyticJacobianMap = { @computeAnalyticDragJacobian, @computeAnalyticWindxJacobian }
%
% This way, the function handles in self.estimatedParamAnalyticJacobianMap are already fully defined
% when calling self.computeAnalyticParamJacobian(), and since self.estimatedParamAnalyticJacobianMap
% is a cell array, it is simply iterated over. Similar reasoning applies for building
% self.estimatedParamSourceMap.
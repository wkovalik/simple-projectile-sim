classdef SensorModel < handle
    properties
        sensorID
        params

        noiseCovar

        projectileModel

        estimatedParams
        estimatedParamJacobianMap
    end

    properties (Dependent)
        invNoiseCovar

        N_PROJECTILE_STATES
        N_ESTIMATED_PARAMS
    end

    % TODO: Try to get this to work. Once it does, move computeParamJacobian() over to here
    properties (Abstract, Constant)
        N_MEASUREMENTS
    end

    methods
        % Constructor ==============================================================================
        function self = SensorModel()
            self.noiseCovar = 0;

            self.estimatedParams = [];
        end

        % Methods ==================================================================================
        function H = computeParamJacobian(self, state)
            H = zeros(self.N_MEASUREMENTS, self.N_ESTIMATED_PARAMS);

            for i = 1:self.N_ESTIMATED_PARAMS
                H(:, i) = self.estimatedParamJacobianMap{i}(self, state);  % See Note 1
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

        function set.estimatedParams(self, estimatedParams)
            if ~isempty(estimatedParams)
                self.estimatedParams = Validator.validateType(estimatedParams, "string");
            end

            self.updateEstimatedParamJacobianMap();
        end

        function updateEstimatedParamJacobianMap(self)
            self.estimatedParamJacobianMap = {};

            if ~isempty(self.estimatedParams)
                for i = 1:self.N_ESTIMATED_PARAMS
                    estimatedParam = self.estimatedParams(i);
    
                    self.estimatedParamJacobianMap{end + 1} = Constants.PARAM_JACOBIAN_MAP(estimatedParam);
                end
            end
        end
    end

    methods (Abstract)
        % Methods ==================================================================================
        y = computeMeasurement(self, state)      
        H = computeStateJacobian(self, state)
    end
end


% Note 1
%
% The original implementation of the computeParamJacobian() function was
%
% | for i = 1:self.N_ESTIMATED_PARAMS
% |     estimatedParam = self.estimatedParams(i);
% | 
% |     A(:, i) = feval(Constants.PARAM_JACOBIAN_MAP(estimatedParam), self, state);
% | end
%
% However, directly looking up the Jacobian function handles in Constants.PARAM_JACOBIAN_MAP ended
% up being quite slow (around twice as slow as the current implementation, in fact). I dunno why.
% My best guess is the local, individual self.computeEstimatedParamJacobian function handles were
% being rebuilt on every single lookup of Constants.PARAM_JACOBIAN_MAP(estimatedParam).
%
% So instead, a local self.estimatedParamJacobianMap property is built which is a cell array of the
% subset of Jacobian function handles in the full Constants.PARAM_JACOBIAN_MAP for the parameters in
% self.estimatedParams only. This cell array lists the handles in the correct order according to the
% order in self.estimatedParams and is rebuilt whenever self.estimatedParams is updated.
%
% For example:
% 
% | self.estimatedParams = ["CD", "vWindx"]
% | self.estimatedParamJacobianMap = { @computeDragJacobian, @computeWindJacobian }
%
% This way, the function handles in self.estimatedParamJacobianMap are already fully defined when
% calling self.computeParamJacobian(), and since self.estimatedParamJacobianMap is a cell array, it
% is simply iterated over.
classdef Propagator < handle
    properties
        propagate
        
        projectileDynamics

        projectile
        earth

        integrator

        stateInterpolant
        stateSTMInterpolant
        paramSTMInterpolant

        interpMethod
        extrapMethod

        estimatedParams

        sensorList
        sensorSampleRecord
    end

    properties (Dependent)
        N_PROJECTILE_STATES
        N_ESTIMATED_PARAMS
        N_SENSORS
    end

    methods
        % Constructor ==============================================================================
        function self = Propagator(projectileDynamics, integrator)
            self.propagate = @(propTime) error("Propagator method has not been selected. Select method by calling Propagator.selectPropagator(method).");

            self.projectileDynamics = projectileDynamics;

            self.projectile = self.projectileDynamics.projectile;
            self.earth = self.projectileDynamics.earth;

            self.integrator = integrator;

            self.interpMethod = Constants.INTERP_METHOD;
            self.extrapMethod = Constants.EXTRAP_METHOD;
            
            self.estimatedParams = [];

            self.sensorList = {};
            self.sensorSampleRecord = dictionary();
        end

        % Propagate methods ========================================================================
        function selectPropagator(self, method)
            switch method
                case "stateOnly"
                    self.propagate = @self.propagateState;
                    self.integrator.stateDerivFn = ...
                        @(varargin) self.projectileDynamics.computeStateDeriv(varargin{:});  % See Note 1
        
                case "stateAndSTM"
                    self.propagate = @self.propagateStateAndSTM;
                    self.integrator.stateDerivFn = ...
                        @(varargin) self.projectileDynamics.computeStateAndSTMDeriv(varargin{:});  % See Note 1
                
                case "stateWithSensors"
                    self.propagate = @self.propagateStateWithSensors;
                    self.integrator.stateDerivFn = ...
                        @(varargin) self.projectileDynamics.computeStateDeriv(varargin{:});  % See Note 1

                otherwise
                    error("Invalid propagation method.")
            end
        end

        function propagateState(self, finalTime)
            timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
            stateHistory = zeros(self.N_PROJECTILE_STATES, Constants.BUFFER_LEN);

            % TODO: Implement a buffer increase function (right now, will throw error if reach history limit)

            time = self.projectile.time;
            state = self.projectile.state.vector;

            timeHistory(1) = time;
            stateHistory(:, 1) = state;
            
            nSteps = 0;
            while time < (finalTime - Constants.TIME_TOLERANCE)  % See Constants:Note 1
                [time, state] = self.integrator.step(time, state);

                nSteps = nSteps + 1;

                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;
            end

            self.projectile.updateState(time, state);

            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];

            self.createStateInterpolant(timeHistory, stateHistory);
        end

        function propagateStateAndSTM(self, finalTime)
            timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
            stateHistory = zeros(self.N_PROJECTILE_STATES, Constants.HISTORY_BUFFER_LEN);
            stateSTMHistory = zeros(self.N_PROJECTILE_STATES ^ 2, Constants.HISTORY_BUFFER_LEN);
            if ~isempty(self.estimatedParams)
                paramSTMHistory = zeros(self.N_PROJECTILE_STATES * self.N_ESTIMATED_PARAMS, Constants.HISTORY_BUFFER_LEN);
            else
                paramSTMHistory = [];
            end

            % TODO: Implement a buffer increase function (right now, will throw error if reach history limit)

            time = self.projectile.time;
            state = self.projectile.state.vector;

            stateSTM = eye(self.N_PROJECTILE_STATES);
            stateSTM = stateSTM(:);
        
            if ~isempty(self.estimatedParams)
                paramSTM = zeros(self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS);
                paramSTM = paramSTM(:);
            else
                paramSTM = [];
            end

            augState = [state; stateSTM; paramSTM];

            timeHistory(1) = time;
            stateHistory(:, 1) = state;
            stateSTMHistory(:, 1) = stateSTM;
            if ~isempty(self.estimatedParams)
                paramSTMHistory(:, 1) = paramSTM;
            end
            
            nSteps = 0;
            while time < (finalTime - Constants.TIME_TOLERANCE)  % See Constants:Note 1
                [time, augState] = self.integrator.step(time, augState);

                state = augState(1:self.N_PROJECTILE_STATES);
                stateSTM = augState(self.N_PROJECTILE_STATES + (1:self.N_PROJECTILE_STATES ^ 2));
                if ~isempty(self.estimatedParams)
                    paramSTM = augState((self.N_PROJECTILE_STATES + self.N_PROJECTILE_STATES ^ 2) + (1:(self.N_PROJECTILE_STATES * self.N_ESTIMATED_PARAMS)));
                end

                nSteps = nSteps + 1;

                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;
                stateSTMHistory(:, nSteps + 1) = stateSTM;
                if ~isempty(self.estimatedParams)
                    paramSTMHistory(:, nSteps + 1) = paramSTM;
                end
            end

            self.projectile.updateState(time, state);

            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];
            stateSTMHistory(:, (nSteps + 2):end) = [];
            if ~isempty(self.estimatedParams)
                paramSTMHistory(:, (nSteps + 2):end) = [];
            end

            self.createStateInterpolant(timeHistory, stateHistory);
            self.createSTMInterpolant(timeHistory, stateSTMHistory, paramSTMHistory);
        end

        function propagateStateWithSensors(self, finalTime)
            timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
            stateHistory = zeros(self.N_PROJECTILE_STATES, Constants.HISTORY_BUFFER_LEN);

            self.initSensorSampleRecord();

            % TODO: Implement a buffer increase function (right now, will throw error if reach history limit)

            time = self.projectile.time;
            state = self.projectile.state.vector;

            self.checkSensors(time, state);

            timeHistory(1) = time;
            stateHistory(:, 1) = state;
            
            nSteps = 0;
            while time < (finalTime - Constants.TIME_TOLERANCE)  % See Constants:Note 1
                [time, state] = self.integrator.step(time, state);

                self.checkSensors(time, state);

                nSteps = nSteps + 1;

                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;
            end

            self.projectile.updateState(time, state);

            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];

            for i = 1:self.N_SENSORS
                sensor = self.sensorList{i};
            
                nSamples = self.sensorSampleRecord(sensor.sensorID).nSamples;
                self.sensorSampleRecord(sensor.sensorID).timeHistory((nSamples + 1):end) = [];
                self.sensorSampleRecord(sensor.sensorID).sampleHistory(:, (nSamples + 1):end) = [];
            end

            self.createStateInterpolant(timeHistory, stateHistory);
        end
        
        % Interpolant methods ======================================================================
        function createStateInterpolant(self, timeHistory, stateHistory)
            self.stateInterpolant = griddedInterpolant(timeHistory, stateHistory', self.interpMethod, self.extrapMethod);
        end

        function state = getState(self, time)
            state = self.stateInterpolant(time)';
        end

        function [timeHistory, stateHistory] = generateStateHistory(self, initTime, finalTime, timeStepOrNumPoints)
            % Mode 1: User specifies time step between points in history
            if timeStepOrNumPoints < (finalTime - initTime)
                timeStep = timeStepOrNumPoints;
                timeHistory = initTime:timeStep:finalTime;

            % Mode 2: User specifies number of points in history
            else
                numOfPoints = timeStepOrNumPoints;
                timeHistory = linspace(initTime, finalTime, numOfPoints);
            end

            stateHistory = self.getState(timeHistory);
        end

        function createSTMInterpolant(self, timeHistory, stateSTMHistory, paramSTMHistory)
            self.stateSTMInterpolant = griddedInterpolant(timeHistory, stateSTMHistory', self.interpMethod, self.extrapMethod);
            if ~isempty(self.estimatedParams)
                self.paramSTMInterpolant = griddedInterpolant(timeHistory, paramSTMHistory', self.interpMethod, self.extrapMethod);
            end
        end

        function [stateSTM, paramSTM] = getSTM(self, time)
            stateSTM = self.stateSTMInterpolant(time)';
            stateSTM = reshape(stateSTM, [self.N_PROJECTILE_STATES, self.N_PROJECTILE_STATES]);

            if ~isempty(self.estimatedParams)
                paramSTM = self.paramSTMInterpolant(time)';
                paramSTM = reshape(paramSTM, [self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS]);
            else
                paramSTM = [];
            end
        end

        % Sensor methods ===========================================================================
        function addSensor(self, sensor)
            self.sensorList(end + 1) = { sensor };
        end

        function checkSensors(self, time, state)
            for i = 1:self.N_SENSORS
                sensor = self.sensorList{i};

                if sensor.checkIfShouldSample(time)
                    sensor.sampleMeasurement(time, state);
                    
                    self.addSampleToRecord(sensor);
                end
            end
        end

        function initSensorSampleRecord(self)
            for i = 1:self.N_SENSORS
                sensor = self.sensorList{i};

                initRecord.nSamples = 0;
                initRecord.timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
                initRecord.sampleHistory = zeros(sensor.N_MEASUREMENTS, Constants.HISTORY_BUFFER_LEN);

                self.sensorSampleRecord = insert(self.sensorSampleRecord, sensor.sensorID, initRecord);
            end
        end

        function addSampleToRecord(self, sensor)
            nSamples = self.sensorSampleRecord(sensor.sensorID).nSamples;
            nSamples = nSamples + 1;

            self.sensorSampleRecord(sensor.sensorID).nSamples = nSamples;
            self.sensorSampleRecord(sensor.sensorID).timeHistory(nSamples) = sensor.sampleTime;
            self.sensorSampleRecord(sensor.sensorID).sampleHistory(:, nSamples) = sensor.sampledMeasurement;
        end

        % Getters ==================================================================================
        function N_PROJECTILE_STATES = get.N_PROJECTILE_STATES(self)
            N_PROJECTILE_STATES = self.projectile.state.N_STATES;
        end

        function N_ESTIMATED_PARAMS = get.N_ESTIMATED_PARAMS(self)
            N_ESTIMATED_PARAMS = length(self.estimatedParams);
        end

        function N_SENSORS = get.N_SENSORS(self)
            N_SENSORS = length(self.sensorList);
        end

        % Setters ==================================================================================
        function set.propagate(self, propagateFn)
            self.propagate = Validator.validateType(propagateFn, "function_handle");
        end

        function set.projectile(self, projectile)
            self.projectile = Validator.validateType(projectile, "Projectile");
        end

        function set.projectileDynamics(self, projectileDynamics)
            self.projectileDynamics = Validator.validateType(projectileDynamics, "ProjectileDynamics");
        end

        function set.integrator(self, integrator)
            self.integrator = Validator.validateType(integrator, "RK4Integrator");
        end

        function set.stateInterpolant(self, stateInterpolant)
            self.stateInterpolant = Validator.validateType(stateInterpolant, "griddedInterpolant");
        end

        function set.interpMethod(self, interpMethod)
            % TODO: Validate from enum of possible values
            self.interpMethod = Validator.validateType(interpMethod, "string");
        end

        function set.extrapMethod(self, extrapMethod)
            % TODO: Validate from enum of possible values
            self.extrapMethod = Validator.validateType(extrapMethod, "string");
        end
    end
end


% Note 1:
%
% |  stateDerivFn = @self.forceModel.computeStateDeriv;
%
% This line should work but doesn't due to how MATLAB handles the @ construction. In particular,
% MATLAB does not create the forceModel.computeStateDeriv() function handle in this scope (where
% self is a Propagator object) but rather in the Integrator.step() scope (where self is the
% Integrator object) since that is where the function is actually invoked. Thus, MATLAB tries
% calling self(Integrator).forceModel.computeStateDeriv(), which obviously does not exist.
%
% |  stateDerivFn = @(varargin) self.forceModel.computeStateDeriv(varargin{:});
%
% This line accomplishes the intended behavior by creating a duplicate function with the exact same
% arguments as forceModel.computeStateDeriv() (using varargin). Thus, the function handle is now
% correctly defined in this scope.
%
% See: https://www.mathworks.com/matlabcentral/answers/489847-how-to-pass-a-method-as-a-function-handle

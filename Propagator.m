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
        sensorSampleHistory
    end

    properties (Dependent)
        N_PROJECTILE_STATES
        N_ESTIMATED_PARAMS
        N_SENSORS
    end

    methods
        % Constructor ==============================================================================
        function self = Propagator(projectileDynamics, integrator)
            % Initialize propagate method (i.e., as unset, so error if try to call without setting)
            self.propagate = @(propTime) error("Propagator method has not been selected. Select method by calling Propagator.selectPropagator(method).");
            
            % Set handle for projectile dynamics
            self.projectileDynamics = projectileDynamics;
            
            % Set handles for projectile and planet
            self.projectile = self.projectileDynamics.projectile;
            self.earth = self.projectileDynamics.earth;

            % Set handle for integrator
            self.integrator = integrator;
            
            % Set initial interpolant methods
            self.interpMethod = Constants.INTERP_METHOD;
            self.extrapMethod = Constants.EXTRAP_METHOD;
            
            % Initialize list of parameters to be estimated (used to compute parameter STM)
            self.estimatedParams = [];
            
            % Initialize list of sensors and measurement histories (used to take sensor measurements)
            self.sensorList = {};
            self.sensorSampleHistory = dictionary();
        end


        % Propagate methods ========================================================================
        function selectPropagator(self, method)
            switch method
                % Propagate projectile state vector only
                case "stateOnly"
                    self.propagate = @self.propagateState;
                    self.integrator.stateDerivFn = ...
                        @(varargin) self.projectileDynamics.computeStateDeriv(varargin{:});  % See Note 1
                
                % Propagate projectile state vector, state STM, and estimated parameter STM (if any)
                case "stateAndSTM"
                    self.propagate = @self.propagateStateAndSTM;
                    self.integrator.stateDerivFn = ...
                        @(varargin) self.projectileDynamics.computeStateAndSTMDeriv(varargin{:});  % See Note 1
                
                % Propagate projectile state vector and take sensor measurements along trajectory
                case "stateWithSensors"
                    self.propagate = @self.propagateStateWithSensors;
                    self.integrator.stateDerivFn = ...
                        @(varargin) self.projectileDynamics.computeStateDeriv(varargin{:});  % See Note 1

                otherwise
                    error("Invalid propagation method.")
            end
        end

        function propagateState(self, finalTime)
            % Preallocate time and projectile state histories
            timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
            stateHistory = zeros(self.N_PROJECTILE_STATES, Constants.BUFFER_LEN);  % TODO: Implement a buffer increase function (right now, will throw error if reach history limit)
            
            % Set initial time and projectile state
            time = self.projectile.time;
            state = self.projectile.state.vector;
            
            % Store initial time and projectile state
            timeHistory(1) = time;
            stateHistory(:, 1) = state;

            % Main propagation loop ----------------------------------------------------------------

            nSteps = 0;
            while time < (finalTime - Constants.TIME_TOLERANCE)  % See Constants:Note 1
                % Step to next time and state
                [nextTime, nextState] = self.integrator.step(time, state);

                % TODO: Implement proper event checking
                % if nextState(3) > 0
                %     break
                % end

                time = nextTime;
                state = nextState;
                
                % Store next time and projectile state
                nSteps = nSteps + 1;
                
                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;
            end

            % --------------------------------------------------------------------------------------
            
            % Update projectile with final time and state
            self.projectile.updateState(time, state);
            
            % Remove all unused entries in time and projectile state histories
            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];
            
            % Create interpolant for projectile trajectory
            self.createStateInterpolant(timeHistory, stateHistory);
        end

        function propagateStateAndSTM(self, finalTime)
            % Preallocate time, state, and STM histories
            timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
            stateHistory = zeros(self.N_PROJECTILE_STATES, Constants.HISTORY_BUFFER_LEN);

            stateSTMHistory = zeros(self.N_PROJECTILE_STATES ^ 2, Constants.HISTORY_BUFFER_LEN);
            if ~isempty(self.estimatedParams)
                paramSTMHistory = zeros(self.N_PROJECTILE_STATES * self.N_ESTIMATED_PARAMS, Constants.HISTORY_BUFFER_LEN);
            else
                paramSTMHistory = [];
            end  % TODO: Implement a buffer increase function (right now, will throw error if reach history limit)
            
            % Set initial time and projectile state
            time = self.projectile.time;
            state = self.projectile.state.vector;
            
            % Set initial state STM
            stateSTM = eye(self.N_PROJECTILE_STATES);
            stateSTM = stateSTM(:);
            
            % Set initial parameter STM (if any estimated parameters)
            if ~isempty(self.estimatedParams)
                paramSTM = zeros(self.N_PROJECTILE_STATES, self.N_ESTIMATED_PARAMS);
                paramSTM = paramSTM(:);
            else
                paramSTM = [];
            end
            
            % Concatenate into augmented initial state
            augState = [state; stateSTM; paramSTM];
            
            % Store initial time, state, and STMs
            timeHistory(1) = time;
            stateHistory(:, 1) = state;

            stateSTMHistory(:, 1) = stateSTM;
            if ~isempty(self.estimatedParams)
                paramSTMHistory(:, 1) = paramSTM;
            end

            % Main propagation loop ----------------------------------------------------------------
            
            nSteps = 0;
            while time < (finalTime - Constants.TIME_TOLERANCE)  % See Constants:Note 1
                % Step to next time and augmented state
                [nextTime, nextAugState] = self.integrator.step(time, augState);
                
                % Deconcatenate augmented state
                nextState = nextAugState(1:self.N_PROJECTILE_STATES);

                nextStateSTM = nextAugState(self.N_PROJECTILE_STATES + (1:self.N_PROJECTILE_STATES ^ 2));
                if ~isempty(self.estimatedParams)
                    nextParamSTM = nextAugState((self.N_PROJECTILE_STATES + self.N_PROJECTILE_STATES ^ 2) + (1:(self.N_PROJECTILE_STATES * self.N_ESTIMATED_PARAMS)));
                end

                % TODO: Implement proper event checking
                % if nextState(3) > 0
                %     break
                % end

                time = nextTime;
                state = nextState;

                stateSTM = nextStateSTM;
                if ~isempty(self.estimatedParams)
                    paramSTM = nextParamSTM;
                end

                augState = nextAugState;

                % Store next time, state, and STMs
                nSteps = nSteps + 1;
                
                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;

                stateSTMHistory(:, nSteps + 1) = stateSTM;
                if ~isempty(self.estimatedParams)
                    paramSTMHistory(:, nSteps + 1) = paramSTM;
                end
            end

            % --------------------------------------------------------------------------------------
            
            % Update projectile with final time and state
            self.projectile.updateState(time, state);
            
            % Remove all unused entries in time, state, and STM histories
            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];

            stateSTMHistory(:, (nSteps + 2):end) = [];
            if ~isempty(self.estimatedParams)
                paramSTMHistory(:, (nSteps + 2):end) = [];
            end
            
            % Create interpolants for projectile trajectory and STM histories
            self.createStateInterpolant(timeHistory, stateHistory);
            self.createSTMInterpolant(timeHistory, stateSTMHistory, paramSTMHistory);
        end

        function propagateStateWithSensors(self, finalTime)
            % Preallocate time and projectile state histories
            timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
            stateHistory = zeros(self.N_PROJECTILE_STATES, Constants.HISTORY_BUFFER_LEN);  % TODO: Implement a buffer increase function (right now, will throw error if reach history limit)

            % Preallocate sensor measurement histories
            self.initSensorSampleHistory();  % TODO: Implement a buffer increase function (right now, will throw error if reach history limit)
            
            % Set initial time and projectile state
            time = self.projectile.time;
            state = self.projectile.state.vector;
            
            % Take initial sensor measurements (if needed)
            self.checkSensors(time, state);
            
            % Store initial time and projectile state
            timeHistory(1) = time;
            stateHistory(:, 1) = state;
            
            % Main propagation loop ----------------------------------------------------------------

            nSteps = 0;
            while time < (finalTime - Constants.TIME_TOLERANCE)  % See Constants:Note 1
                % Step to next time and state
                [nextTime, nextState] = self.integrator.step(time, state);

                % TODO: Implement proper event checking
                if nextState(3) > 0
                    break
                end

                time = nextTime;
                state = nextState;
                
                % Take sensor measurements (if needed)
                self.checkSensors(time, state);
                
                % Store next time and projectile state
                nSteps = nSteps + 1;

                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;
            end

            % --------------------------------------------------------------------------------------
            
            % Update projectile with final time and state
            self.projectile.updateState(time, state);
            
            % Remove all unused entries in time and projectile state histories
            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];

            % Remove all unused entries in sensor measurement histories
            for i = 1:self.N_SENSORS
                sensor = self.sensorList{i};
            
                nSamples = self.sensorSampleHistory(sensor.sensorID).nSamples;
                self.sensorSampleHistory(sensor.sensorID).timeHistory((nSamples + 1):end) = [];
                self.sensorSampleHistory(sensor.sensorID).sampleHistory(:, (nSamples + 1):end) = [];
            end
            
            % Create interpolant for projectile trajectory
            self.createStateInterpolant(timeHistory, stateHistory);
        end
        

        % Interpolant methods ======================================================================
        function createStateInterpolant(self, timeHistory, stateHistory)
            % Creates interpolant from state history array

            self.stateInterpolant = griddedInterpolant(timeHistory, stateHistory', self.interpMethod, self.extrapMethod);
        end

        function state = getState(self, time)
            % Gets interpolated state at desired time

            state = self.stateInterpolant(time)';
        end

        function [timeHistory, stateHistory] = generateStateHistory(self, initTime, finalTime, timeStepOrNumPoints)
            % Creates interpolated state history array

            % Mode 1: User specifies time step between points in history
            if timeStepOrNumPoints < (finalTime - initTime)
                timeStep = timeStepOrNumPoints;
                timeHistory = initTime:timeStep:finalTime;

            % Mode 2: User specifies number of points in history (useful for plotting)
            else
                numOfPoints = timeStepOrNumPoints;
                timeHistory = linspace(initTime, finalTime, numOfPoints);
            end

            stateHistory = self.getState(timeHistory);
        end

        function createSTMInterpolant(self, timeHistory, stateSTMHistory, paramSTMHistory)
            % Creates interpolants of vectorized STMs from vectorized STM history arrays

            self.stateSTMInterpolant = griddedInterpolant(timeHistory, stateSTMHistory', self.interpMethod, self.extrapMethod);
            if ~isempty(self.estimatedParams)
                self.paramSTMInterpolant = griddedInterpolant(timeHistory, paramSTMHistory', self.interpMethod, self.extrapMethod);
            end
        end

        function [stateSTM, paramSTM] = getSTM(self, time)
            % Gets interpolated matrix STMs at desired times

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
            % Add new sensor to sensor list

            self.sensorList(end + 1) = { sensor };
        end

        function checkSensors(self, time, state)
            % Determines if any sensors are ready to take a new measurement, and if so, takes those measurements

            for i = 1:self.N_SENSORS
                sensor = self.sensorList{i};

                if sensor.checkIfShouldSample(time)
                    % Take new sensor measurement
                    sensor.sampleMeasurement(time, state);
                    
                    % Add sensor measurement to respective history (according to sensor ID)
                    self.addSampleToHistory(sensor);
                end
            end
        end

        function initSensorSampleHistory(self)
            % Preallocate histories to store respective sensor measurements (according to sensor ID)

            for i = 1:self.N_SENSORS
                sensor = self.sensorList{i};

                initHistory.nSamples = 0;
                initHistory.timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
                initHistory.sampleHistory = zeros(sensor.N_MEASUREMENTS, Constants.HISTORY_BUFFER_LEN);

                self.sensorSampleHistory = insert(self.sensorSampleHistory, sensor.sensorID, initHistory);
            end
        end

        function addSampleToHistory(self, sensor)
            % Add new sensor measurement to respective history (according to sensor ID)

            nSamples = self.sensorSampleHistory(sensor.sensorID).nSamples;
            nSamples = nSamples + 1;

            self.sensorSampleHistory(sensor.sensorID).nSamples = nSamples;
            self.sensorSampleHistory(sensor.sensorID).timeHistory(nSamples) = sensor.sampleTime;
            self.sensorSampleHistory(sensor.sensorID).sampleHistory(:, nSamples) = sensor.sampledMeasurement;
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

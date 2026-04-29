classdef Propagator < handle
    % TODO: Implement a history length increase function (right now, will throw error if reach history limit)
    % TODO: Implement proper event checking

    properties
        projectileDynamics
        integrator
    end


    methods
        % Constructor ==============================================================================

        function self = Propagator(projectileDynamics, integrator)
            % Set handles for projectile dynamics and integrator
            if nargin == 2
                self.projectileDynamics = projectileDynamics;
                self.integrator = integrator;

            elseif nargin == 1
                self.projectileDynamics = projectileDynamics;
                self.integrator = Integrator();

            else
                error("Not enough input arguments. Requires at least projectileDynamics.")
            end
        end


        % Propagate methods ========================================================================
        
        function [timeHistory, stateHistory] = propagate(self, finalTime)
            self.projectileDynamics.update();

            % Set integrator derivative function
            self.integrator.computeStateDeriv = @(varargin) self.projectileDynamics.computeStateDeriv(varargin{:});  % See Note 1

            % Preallocate time and state histories
            nHistoryLen = Settings.DEFAULT_HISTORY_LEN;
            nStates = self.projectileDynamics.projectile.nStates;

            timeHistory = zeros(1, nHistoryLen);
            stateHistory = zeros(nStates, nHistoryLen);
            
            % Set initial time and state
            time = self.projectileDynamics.projectile.time;
            state = self.projectileDynamics.projectile.state;
            
            % Store initial time and state
            timeHistory(1) = time;
            stateHistory(:, 1) = state;
            
            % --------------------------------------------------------------------------------------
            % Begin propagation loop
            % --------------------------------------------------------------------------------------

            nSteps = 0;
            while time < (finalTime - Settings.DEFAULT_TIME_TOL)  % See Settings: Note 1
                % Step to next time and state
                [nextTime, nextState] = self.integrator.step(time, state);

                % Detect ground impact
                if (nextTime > 5) && (nextState(3) > 0)
                    break
                end

                time = nextTime;
                state = nextState;
                
                % Store next time and state
                nSteps = nSteps + 1;
                
                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;
            end

            % --------------------------------------------------------------------------------------
            % End propagation loop
            % --------------------------------------------------------------------------------------
            
            % Remove all unused entries in time and state histories
            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];
        end


        function [timeHistory, stateHistory, stateSTMHistory, paramSTMHistory] = propagateWithSTM(self, finalTime)
            self.projectileDynamics.update();

            % Set integrator derivative function
            self.integrator.computeStateDeriv = @(varargin) self.projectileDynamics.computeAugStateDeriv(varargin{:});  % See Note 1
            
            % Preallocate time and state histories
            nHistoryLen = Settings.DEFAULT_HISTORY_LEN;
            nStates = self.projectileDynamics.projectile.nStates;
            nEstimatedParams = self.projectileDynamics.projectile.nEstimatedParams + self.projectileDynamics.planet.nEstimatedParams;

            includeParamSTM = logical(nEstimatedParams);

            if ~includeParamSTM
                paramSTM = [];
                paramSTMHistory = [];
            end

            timeHistory = zeros(1, nHistoryLen);
            stateHistory = zeros(nStates, nHistoryLen);
            
            % Preallocate STM histories
            stateSTMHistory = eye(nStates ^ 2, nHistoryLen);
            if includeParamSTM
                paramSTMHistory = eye(nStates * nEstimatedParams, nHistoryLen);
            end

            % Set initial time and state
            time = self.projectileDynamics.projectile.time;
            state = self.projectileDynamics.projectile.state;
            
            % Set initial STMs
            stateSTM = eye(nStates);
            stateSTM = stateSTM(:);
            
            if includeParamSTM
                paramSTM = zeros(nStates, nEstimatedParams);
                paramSTM = paramSTM(:);
            end
            
            % Create augmented state
            augState = [state; stateSTM; paramSTM];
            
            % Store initial time and state
            timeHistory(1) = time;
            stateHistory(:, 1) = state;

            % Store initial STMs
            stateSTMHistory(:, 1) = stateSTM;
            if includeParamSTM
                paramSTMHistory(:, 1) = paramSTM;
            end

            % --------------------------------------------------------------------------------------
            % Begin propagation loop
            % --------------------------------------------------------------------------------------

            iStateEnd = nStates;
            iStateSTMEnd = nStates + nStates ^ 2;

            nSteps = 0;
            while time < (finalTime - Settings.DEFAULT_TIME_TOL)  % See Settings: Note 1
                % Step to next time and augmented state
                [nextTime, nextAugState] = self.integrator.step(time, augState);
                
                % Extract next state
                nextState = nextAugState(1:iStateEnd);

                % Detect ground impact
                if (nextTime > 5) && (nextState(3) > 0)
                    break
                end

                time = nextTime;
                state = nextState;
                augState = nextAugState;

                % Extract STMs
                stateSTM = augState((iStateEnd + 1):iStateSTMEnd);
                if includeParamSTM
                    paramSTM = augState((iStateSTMEnd + 1):end);
                end
                
                % Store next time and state
                nSteps = nSteps + 1;
                
                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;

                % Store STMs
                stateSTMHistory(:, nSteps + 1) = stateSTM;
                if includeParamSTM
                    paramSTMHistory(:, nSteps + 1) = paramSTM;
                end
            end

            % --------------------------------------------------------------------------------------
            % End propagation loop
            % --------------------------------------------------------------------------------------
            
            % Remove all unused entries in time and state histories
            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];

            % Remove all unused entries in STM histories
            stateSTMHistory(:, (nSteps + 2):end) = [];
            if includeParamSTM
                paramSTMHistory(:, (nSteps + 2):end) = [];
            end
        end


        function [timeHistory, stateHistory, measHistory] = propagateWithSensors(self, finalTime, sensorArray)
            self.projectileDynamics.update();
            
            % Throw warning if integrator step size is too coarse for sensor sampling rates
            if Utils.isIntegratorStepTooCoarse(sensorArray, self.integrator.stepPeriod)
                warning("Integrator step size is too coarse. At least one sensor will not accurately sample on time.")
            end
            
            % Set integrator derivative function
            self.integrator.computeStateDeriv = @(varargin) self.projectileDynamics.computeStateDeriv(varargin{:});  % See Note 1
            
            % Preallocate time and state histories
            nHistoryLen = Settings.DEFAULT_HISTORY_LEN;
            nStates = self.projectileDynamics.projectile.nStates;

            timeHistory = zeros(1, nHistoryLen);
            stateHistory = zeros(nStates, nHistoryLen);
            
            % Set initial time and state
            time = self.projectileDynamics.projectile.time;
            state = self.projectileDynamics.projectile.state;
            
            % Store initial time and state
            timeHistory(1) = time;
            stateHistory(:, 1) = state;

            % Preallocate sensor measurement history
            nMeasMax = Utils.findLongestSensorMeas(sensorArray);
            measHistory = zeros(2 + nMeasMax, nHistoryLen);
            
            % Take any initial sensor measurement(s)
            nSamples = 0;
            for i = 1:length(sensorArray)
                sensor = sensorArray{i};

                if sensor.shouldTakeMeasurement(time)
                    sample = sensor.takeMeasurement(state);
                    
                    % Store sensor measurement(s)
                    nSamples = nSamples + 1;

                    nMeas = sensor.nMeas;
                    iMeasEnd = 3 + (nMeas - 1);

                    measHistory(1, nSamples) = sensor.ID;
                    measHistory(2, nSamples) = time;
                    measHistory(3:iMeasEnd, nSamples) = sample;
                end
            end

            % --------------------------------------------------------------------------------------
            % Begin propagation loop
            % --------------------------------------------------------------------------------------

            nSteps = 0;
            while time < (finalTime - Settings.DEFAULT_TIME_TOL)  % See Settings: Note 1
                % Step to next time and state
                [nextTime, nextState] = self.integrator.step(time, state);

                % Detect ground impact
                if (nextTime > 5) && (nextState(3) > 0)
                    break
                end

                time = nextTime;
                state = nextState;
                
                % Take any sensor measurement(s)
                for i = 1:length(sensorArray)
                    sensor = sensorArray{i};

                    if sensor.shouldTakeMeasurement(time)
                        sample = sensor.takeMeasurement(state);
                        
                        % Store sensor measurement(s)
                        nSamples = nSamples + 1;

                        nMeas = sensor.nMeas;
                        iMeasEnd = 3 + (nMeas - 1);

                        measHistory(1, nSamples) = sensor.ID;
                        measHistory(2, nSamples) = time;
                        measHistory(3:iMeasEnd, nSamples) = sample;
                    end
                end
                
                % Store next time and state
                nSteps = nSteps + 1;
                
                timeHistory(nSteps + 1) = time;
                stateHistory(:, nSteps + 1) = state;
            end

            % --------------------------------------------------------------------------------------
            % End propagation loop
            % --------------------------------------------------------------------------------------
            
            % Remove all unused entries in time and state histories
            timeHistory((nSteps + 2):end) = [];
            stateHistory(:, (nSteps + 2):end) = [];

            % Remove all unused entries in sensor measurement history
            measHistory(:, (nSamples + 1):end) = [];
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

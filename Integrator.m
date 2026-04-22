classdef Integrator < handle
    properties
        stepPeriod

        computeStateDeriv
    end
    
    properties (Constant)
        DEFAULT_STEP_PERIOD = 0.01;  % (s)
    end

    methods
        % Constructor ==============================================================================

        function self = Integrator(stepPeriod)
            % Set integration step size
            if nargin == 1
                self.stepPeriod = stepPeriod;
            else
                self.stepPeriod = self.DEFAULT_STEP_PERIOD;
            end
        end
        

        % Methods ==================================================================================

        function [time, state] = step(self, time, state)
            % Compute intermediate state derivatives
            k1 = self.computeStateDeriv(state);
            k2 = self.computeStateDeriv(state + k1 * self.stepPeriod / 2);
            k3 = self.computeStateDeriv(state + k2 * self.stepPeriod / 2);
            k4 = self.computeStateDeriv(state + k3 * self.stepPeriod);
            
            % Compute next state and time
            state = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6 * self.stepPeriod;
            time = time + self.stepPeriod;
        end
        

        % Setters ==================================================================================
        
        function set.stepPeriod(self, stepPeriod)
            if Constants.VALIDATE_FLAG
                self.stepPeriod = Validator.validateType(stepPeriod, "double");
            else
                self.stepPeriod = stepPeriod;
            end
        end

        function set.computeStateDeriv(self, stateDerivFn)
            if Constants.VALIDATE_FLAG
                self.computeStateDeriv = Validator.validateType(stateDerivFn, "function_handle");
            else
                self.computeStateDeriv = stateDerivFn;
            end
        end
    end
end
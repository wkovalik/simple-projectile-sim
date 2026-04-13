classdef RK4Integrator < handle
    properties
        stepPeriod

        stateDerivFn
    end
    
    properties (Constant)
        DEFAULT_STEP_PERIOD = 0.005;
    end

    methods
        % Constructor ==============================================================================
        function self = RK4Integrator(initStepPeriod)
            self.stepPeriod = self.DEFAULT_STEP_PERIOD;

            if nargin > 0
                self.stepPeriod = initStepPeriod;
            end
        end
        
        % Methods ==================================================================================
        function [time, state] = step(self, time, state)
            k1 = self.stateDerivFn(state);
            k2 = self.stateDerivFn(state + k1 * self.stepPeriod / 2);
            k3 = self.stateDerivFn(state + k2 * self.stepPeriod / 2);
            k4 = self.stateDerivFn(state + k3 * self.stepPeriod);

            state = state + (k1 + 2 * k2 + 2 * k3 + k4) / 6 * self.stepPeriod;
            time = time + self.stepPeriod;
        end
        
        % Setters ==================================================================================
        function set.stepPeriod(self, stepPeriod)
            self.stepPeriod = Validator.validateType(stepPeriod, "double");
        end

        function set.stateDerivFn(self, stateDerivFn)
            self.stateDerivFn = Validator.validateType(stateDerivFn, "function_handle");
        end
    end
end
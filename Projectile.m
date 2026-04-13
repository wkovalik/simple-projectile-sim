classdef Projectile < handle
    properties        
        time
        state

        params
    end

    methods
        % Constructor ==============================================================================
        function self = Projectile()            
            self.time = 0;  % (s)
            self.state = ProjectileState();

            self.params.m = Param();   % Mass             (kg)
            self.params.S = Param();   % Reference area   (m^2)
            self.params.CD = Param();  % Drag coefficient
        end

        % Methods ==================================================================================
        function updateState(self, time, stateVector, stateCovar)
            self.time = time;
            self.state.vector = stateVector;

            if nargin > 3
                self.state.covar = stateCovar;
            end
        end

        % Setters ==================================================================================
        function set.time(self, time)
            self.time = Validator.validateType(time, "double");
        end

        function set.state(self, state)
            self.state = Validator.validateType(state, "ProjectileState");
        end

        function set.params(self, params)
            self.params = Validator.validateStructArrayTypes(params, "Param");
        end
    end
end
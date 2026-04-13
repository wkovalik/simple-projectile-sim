classdef ProjectileState < handle
    properties
        vector
        covar
    end

    properties (Dependent)
        x   % Position x-component (m)
        y   % Position y-component (m)
        z   % Position z-component (m)
        vx  % Velocity x-component (m/s)
        vy  % Velocity y-component (m/s)
        vz  % Velocity z-component (m/s)

        invCovar
    end

    properties (Constant)
        N_STATES = 6;
    end

    methods
        % Constructor ==============================================================================
        function self = ProjectileState(initVector, initCovar)
            self.vector = zeros(self.N_STATES, 1);
            self.covar = zeros(self.N_STATES);

            if nargin > 0
                self.vector = initVector;
            end
            
            if nargin > 1
                self.covar = initCovar;
            end
        end
        
        % Getters ==================================================================================
        function x = get.x(self)
            x.value = self.vector(1);
            x.covar = self.covar(1, 1);
        end

        function y = get.y(self)
            y.value = self.vector(2);
            y.covar = self.covar(2, 2);
        end

        function z = get.z(self)
            z.value = self.vector(3);
            z.covar = self.covar(3, 3);
        end

        function vx = get.vx(self)
            vx.value = self.vector(4);
            vx.covar = self.covar(4, 4);
        end

        function vy = get.vy(self)
            vy.value = self.vector(5);
            vy.covar = self.covar(5, 5);
        end

        function vz = get.vz(self)
            vz.value = self.vector(6);
            vz.covar = self.covar(6, 6);
        end

        function invCovar = get.invCovar(self)
            invCovar = inv(self.covar);
        end

        % Setters ==================================================================================
        function set.vector(self, vector)
            vector = Validator.validateType(vector, "double");
            self.vector = Validator.validateSize(vector, [self.N_STATES, 1]);
        end

        function set.covar(self, covar)
            covar = Validator.validateType(covar, "double");
            self.covar = Validator.validateSize(covar, [self.N_STATES, self.N_STATES]);
        end
    end
end
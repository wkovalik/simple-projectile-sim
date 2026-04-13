classdef ProjectileState < handle
    properties
        vector
        covar
    end

    properties (Dependent)
        x
        y
        vx
        vy

        invCovar
    end

    properties (Constant)
        N_STATES = 4;
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

        function vx = get.vx(self)
            vx.value = self.vector(3);
            vx.covar = self.covar(3, 3);
        end

        function vy = get.vy(self)
            vy.value = self.vector(4);
            vy.covar = self.covar(4, 4);
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
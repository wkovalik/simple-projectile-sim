classdef Earth < handle
    properties
        params
    end

    properties (Constant)
        DEFAULT_G = 9.81;      % (m/s^2)
        DEFAULT_RHO0 = 1.225;  % (kg/m^3)
        DEFAULT_H = 8500;      % (m)
    end

    methods
        % Constructor ==============================================================================
        function self = Earth()
            self.params.g = Param(self.DEFAULT_G);      % Gravitational acceleration (m/s^2)

            self.params.rho0 = Param(self.DEFAULT_RHO0);  % Atmospheric density at surface (kg/m^3)
            self.params.H = Param(self.DEFAULT_H);        % Atmospheric scale height (m)

            self.params.vWindx = Param();               % Wind velocity x-component (m/s)
            self.params.vWindy = Param();               % Wind velocity y-component (m/s)
        end


        % Methods ==================================================================================
        function rho = computeAtmosphericDensity(self, h)
            rho0 = self.params.rho0.value;
            H = self.params.H.value;

            rho = rho0 * exp(-h / H);
        end


        % Setters ==================================================================================
        function set.params(self, params)
            self.params = Validator.validateStructArrayTypes(params, "Param");
        end
    end
end
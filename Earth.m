classdef Earth < handle
    properties
        params
    end

    properties (Constant)
        DEFAULT_G = 9.81;     % (m/s^2)
        DEFAULT_RHO = 1.225;  % (kg/m^3)
    end

    methods
        % Constructor ==============================================================================
        function self = Earth()
            self.params.g = Param(self.DEFAULT_G);      % Gravitational acceleration (m/s^2)
            self.params.rho = Param(self.DEFAULT_RHO);  % Atmospheric density (kg/m^3)
            self.params.vWindx = Param();               % Wind velocity x-component (m/s)
            self.params.vWindy = Param();               % Wind velocity y-component (m/s)
        end


        % Setters ==================================================================================
        function set.params(self, params)
            self.params = Validator.validateStructArrayTypes(params, "Param");
        end
    end
end
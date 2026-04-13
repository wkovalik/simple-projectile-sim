classdef Earth < handle
    properties
        params
    end

    properties (Constant)
        DEFAULT_G = 9.81;
        DEFAULT_RHO = 1.225;
    end

    methods
        % Constructor ==============================================================================
        function self = Earth()
            self.params.g = Param(self.DEFAULT_G);
            self.params.rho = Param(self.DEFAULT_RHO);
            self.params.vWindx = Param();
        end

        % Setters ==================================================================================
        function set.params(self, params)
            self.params = Validator.validateStructArrayTypes(params, "Param");
        end
    end
end
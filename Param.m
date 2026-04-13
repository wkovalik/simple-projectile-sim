classdef Param < handle
    properties
        value
        covar

        isEstimated
    end

    methods
        % Constructor ==============================================================================
        function self = Param(value, covar)
            self.value = 0;
            self.covar = 0;

            if nargin > 0
                self.value = value;
            end
            if nargin > 1
                self.covar = covar;
            end

            self.isEstimated = false;
        end


        % Setters ==================================================================================
        function set.value(self, value)
            self.value = Validator.validateType(value, "double");
        end

        function set.covar(self, covar)
            self.covar = Validator.validateType(covar, "double");
        end

        function set.isEstimated(self, isEstimatedFlag)
            self.isEstimated = Validator.validateType(isEstimatedFlag, "logical");
        end
    end
end
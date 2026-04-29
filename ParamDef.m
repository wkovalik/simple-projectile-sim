classdef ParamDef
    properties
        value = 0;
        covar = 0;
        isEstimated = false;
    end


    methods
        % Constructor ==============================================================================

        function obj = ParamDef(value, covar, isEstimated)
            if nargin == 3
                obj.value = value;
                obj.covar = covar;
                obj.isEstimated = isEstimated;

            elseif nargin == 2
                obj.value = value;
                obj.covar = covar;

            elseif nargin == 1
                obj.value = value;

            end
        end

        
        % Setters ==================================================================================
        
        function obj = set.value(obj, value)
            if Settings.VALIDATE_FLAG
                obj.value = Validator.validateType(value, "double");
            else
                obj.value = value;
            end
        end

        function obj = set.covar(obj, covar)
            if Settings.VALIDATE_FLAG
                obj.covar = Validator.validateType(covar, "double");
            else
                obj.covar = covar;
            end
        end

        function obj = set.isEstimated(obj, isEstimated)
            if Settings.VALIDATE_FLAG
                obj.isEstimated = Validator.validateType(isEstimated, "logical");
            else
                obj.isEstimated = isEstimated;
            end
        end
    end
end
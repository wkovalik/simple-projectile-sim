classdef ParamTableDef
    properties
        xValues = 0;
        yValues = 0;
        yCovars = 0;
        yIsEstimated = false;
    end

    properties (Dependent)
        nValues
    end


    methods
        % Constructor ==============================================================================

        function obj = ParamTableDef(xValues, yValues, yCovars, yIsEstimated)
            if nargin == 4
                obj.xValues = xValues;
                obj.yValues = yValues;
                obj.yCovars = yCovars;
                obj.yIsEstimated = yIsEstimated;
            
            elseif nargin == 3
                obj.xValues = xValues;
                obj.yValues = yValues;
                obj.yCovars = yCovars;
                obj.yIsEstimated = false(obj.nValues, 1);

            elseif nargin == 2
                obj.xValues = xValues;
                obj.yValues = yValues;
                obj.yCovars = zeros(obj.nValues, 1);
                obj.yIsEstimated = false(obj.nValues, 1);

            elseif nargin == 1
                obj.xValues = xValues;
                obj.yValues = zeros(obj.nValues, 1);
                obj.yCovars = zeros(obj.nValues, 1);
                obj.yIsEstimated = false(obj.nValues, 1);
            
            else
                error("Not enough input parameters.")
            end
        end

        
        % Getters ==================================================================================

        function nValues = get.nValues(obj)
            nValues = length(obj.xValues);
        end

        
        % Setters ==================================================================================
        
        function obj = set.xValues(obj, xValues)
            if Constants.VALIDATE_FLAG
                obj.xValues = Validator.validateType(xValues, "double");
            else
                obj.xValues = xValues;
            end
            
            % Reset all other fields if table length is changed
            if obj.nValues ~= length(obj.yValues)
                obj.yValues = zeros(obj.nValues, 1);
                obj.yCovars = zeros(obj.nValues, 1);
                obj.yIsEstimated = false(obj.nValues, 1);
            end
        end

        function obj = set.yValues(obj, yValues)
            if Constants.VALIDATE_FLAG
                yValues = Validator.validateType(yValues, "double");
                obj.yValues = Validator.validateSize(yValues, [obj.nValues, 1]);
            else
                obj.yValues = yValues;
            end
        end

        function obj = set.yCovars(obj, yCovars)
            if Constants.VALIDATE_FLAG
                yCovars = Validator.validateType(yCovars, "double");
                obj.yCovars = Validator.validateSize(yCovars, [obj.nValues, 1]);
            else
                obj.yCovars = yCovars;
            end
        end

        function obj = set.yIsEstimated(obj, yIsEstimated)
            if Constants.VALIDATE_FLAG
                yIsEstimated = Validator.validateType(yIsEstimated, "logical");
                obj.yIsEstimated = Validator.validateSize(yIsEstimated, [obj.nValues, 1]);
            else
                obj.yIsEstimated = yIsEstimated;
            end
        end
    end
end
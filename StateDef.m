classdef StateDef
    properties
        time = 0;
        state = [];
        covar = [];
    end

    properties (Constant)
        nStates = 6;
    end


    methods
        % Constructor ==============================================================================

        function obj = StateDef(time, state, covar)
            if nargin == 3
                obj.time = time;
                obj.state = state;
                obj.covar = covar;

            elseif nargin == 2
                obj.time = time;
                obj.state = state;
                obj.covar = zeros(obj.nStates);

            elseif nargin == 1
                obj.time = time;
                obj.state = zeros(obj.nStates, 1);
                obj.covar = zeros(obj.nStates);
            
            elseif nargin == 0
                obj.time = 0;
                obj.state = zeros(obj.nStates, 1);
                obj.covar = zeros(obj.nStates);

            end
        end

        
        % Setters ==================================================================================

        function obj = set.time(obj, time)
            if Settings.VALIDATE_FLAG
                obj.time = Validator.validateType(time, "double");
            else
                obj.time = time;
            end
        end

        function obj = set.state(obj, state)
            if Settings.VALIDATE_FLAG
                state = Validator.validateType(state, "double");
                obj.state = Validator.validateSize(state, [obj.nStates, 1]);
            else
                obj.state = state;
            end
        end

        function obj = set.covar(obj, covar)
            if Settings.VALIDATE_FLAG
                covar = Validator.validateType(covar, "double");
                obj.covar = Validator.validateSize(covar, [obj.nStates, obj.nStates]);
            else
                obj.covar = covar;
            end
        end
    end
end
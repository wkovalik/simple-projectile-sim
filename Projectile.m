classdef Projectile < handle
    % TODO: Make a static Kernel class. Then set kernel property to choose which kernel to use

    properties
        stateDef
        paramDefs
        
        time = 0;
        state = [];
        stateCovar = [];

        params = [];

        estimatedParams = [];
        estimatedParamCovar = [];

        aeroModel
    end

    properties (SetAccess = private)
        estimatedParamIdxs = [];

        mIdx = 0;
        SIdx = 0;

        CDIdx = 0;
        CDTable_Mach0Idx = 0;
        CDTable_CD0Idx = 0;
        CDTable_Len = 0;

        computeAeroCoeffs
    end

    properties (Dependent, SetAccess = private)
        nParams
        nEstimatedParams
    end

    properties (Constant)
        nStates = 6;

        xIdx = 1;
        yIdx = 2;
        zIdx = 3;
        vxIdx = 4;
        vyIdx = 5;
        vzIdx = 6;

        DEFAULT_M = 0.145;
        % DEFAULT_D = 0.075;
        DEFAULT_S = 0.004417865;  % TODO: Compute using d as param instead

        DEFAULT_CD = 0.15;
        DEFAULT_CD_TABLE_X = [0; 1];
        DEFAULT_CD_TABLE_Y = [0.15; 0.15];

        VALID_AERO_MODELS = ["constant", "table"];
    end
    

    methods
        % Constructor ==============================================================================

        function self = Projectile(aeroModel)
            self.stateDef = StateDef();

            self.updateState();

            self.paramDefs.m = ParamDef(self.DEFAULT_M);
            self.paramDefs.S = ParamDef(self.DEFAULT_S);
            
            if nargin == 0
                self.aeroModel = "constant";
            elseif nargin == 1
                self.aeroModel = aeroModel;
            else
                error("Too many input parameters.")
            end

            self.updateModels();
        end

        
        % Update methods ===========================================================================

        function updateState(self)
            self.time = self.stateDef.time;
            self.state = self.stateDef.state;
            self.stateCovar = self.stateDef.covar;
        end


        function updateModels(self)
            switch self.aeroModel
                case "constant"
                    self.computeAeroCoeffs = @self.constantAeroModel;

                    self.paramDefs.CD = ParamDef(self.DEFAULT_CD);

                case "table"
                    self.computeAeroCoeffs = @self.tableAeroModel;

                    self.paramDefs.CD = ParamTableDef(self.DEFAULT_CD_TABLE_X, self.DEFAULT_CD_TABLE_Y);

                otherwise
                    error("Invalid aerodynamics model: %s.", self.aeroModel)
            end

            self.updateParams();
        end


        function updateParams(self)
            % See Note 1 regarding parameter indexing

            self.params = [];
            
            self.mIdx = 0;
            self.SIdx = 0;
    
            self.CDIdx = 0;
            self.CDTable_Mach0Idx = 0;
            self.CDTable_CD0Idx = 0;
            self.CDTable_Len = 0;
            
            % m
            self.mIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.m.value];
            
            % S
            self.SIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.S.value];
            
            switch self.aeroModel
                case "constant"
                    % CD
                    self.CDIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CD.value];

                case "table"
                    % CD
                    self.CDTable_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CD.xValues];

                    self.CDTable_CD0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CD.yValues];
                    
                    self.CDTable_Len = self.paramDefs.CD.nValues;
            end

            self.updateEstimatedParams();
        end


        function updateEstimatedParams(self)
            self.estimatedParams = [];
            self.estimatedParamCovar = [];
            self.estimatedParamIdxs = [];

            % m
            if self.paramDefs.m.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.m.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.m.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.mIdx];
            end
            
            % S
            if self.paramDefs.S.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.S.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.S.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.SIdx];
            end
            
            switch self.aeroModel
                case "constant"
                    % CD
                    if self.paramDefs.CD.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CD.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CD.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CDIdx];
                    end

                case "table"
                    % CD
                    for i = 1:length(self.paramDefs.CD.yValues)
                        if self.paramDefs.CD.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CD.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CD.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CDTable_CD0Idx + (i - 1)];
                        end
                    end
            end
        end

        
        % Model methods ============================================================================

        function CD = constantAeroModel(self, ~)
            CD = self.params(self.CDIdx);
        end


        function CD = tableAeroModel(self, mach)
            dmach = self.params(self.CDTable_Mach0Idx + 1) - self.params(self.CDTable_Mach0Idx);

            CD = 0;
            for i = 0:(self.CDTable_Len - 1)
                mach_i = self.params(self.CDTable_Mach0Idx + i);
                CD_i = self.params(self.CDTable_CD0Idx + i);

                CD = CD + CD_i * self.linearKernel(mach - mach_i, dmach);
            end
        end
        

        function k = linearKernel(~, x, dx)
            xx = x / dx;

            if (-1 < xx) && (xx < 0)
                k = 1 + xx;
            elseif (0 <= xx) && (xx < 1)
                k = 1 - xx;
            else
                k = 0;
            end
        end

        
        % Getters ==================================================================================

        function nParams = get.nParams(self)
            nParams = length(self.params);
        end

        function nEstimatedParams = get.nEstimatedParams(self)
            nEstimatedParams = length(self.estimatedParams);
        end

        
        % Setters ==================================================================================

        function set.stateDef(self, stateDef)
            if Constants.VALIDATE_FLAG
                self.stateDef = Validator.validateType(stateDef, "StateDef");
            else
                self.stateDef = stateDef;
            end
        end

        function set.paramDefs(self, paramDefs)
            if Constants.VALIDATE_FLAG
                self.paramDefs = Validator.validateFieldTypes(paramDefs, ["ParamDef", "ParamTableDef"]);
            else
                self.paramDefs = paramDefs;
            end
        end

        function set.state(self, state)
            if Constants.VALIDATE_FLAG
                state = Validator.validateType(state, "double");
                self.state = Validator.validateSize(state, [self.nStates, 1]);
            else
                self.state = state;
            end
        end

        function set.stateCovar(self, stateCovar)
            if Constants.VALIDATE_FLAG
                stateCovar = Validator.validateType(stateCovar, "double");
                self.stateCovar = Validator.validateSize(stateCovar, [self.nStates, self.nStates]);
            else
                self.stateCovar = stateCovar;
            end
        end

        function set.params(self, params)
            if Constants.VALIDATE_FLAG
                self.params = Validator.validateType(params, "double");
            else
                self.params = params;
            end
        end

        function set.estimatedParams(self, estimatedParams)
            if Constants.VALIDATE_FLAG
                self.estimatedParams = Validator.validateType(estimatedParams, "double");
            else
                self.estimatedParams = estimatedParams;
            end
        end

        function set.estimatedParamCovar(self, estimatedParamCovar)
            if Constants.VALIDATE_FLAG
                estimatedParamCovar = Validator.validateType(estimatedParamCovar, "double");
                self.estimatedParamCovar = Validator.validateSize(estimatedParamCovar, [self.nEstimatedParams, self.nEstimatedParams]);
            else
                self.estimatedParamCovar = estimatedParamCovar;
            end
        end

        function set.estimatedParamIdxs(self, estimatedParamIdxs)
            if Constants.VALIDATE_FLAG
                self.estimatedParamIdxs = Validator.validateType(estimatedParamIdxs, "double");
            else
                self.estimatedParamIdxs = estimatedParamIdxs;
            end
        end

        function set.aeroModel(self, aeroModel)
            if Constants.VALIDATE_FLAG
                self.aeroModel = Validator.validateString(aeroModel, self.VALID_AERO_MODELS);
            else
                self.aeroModel = aeroModel;
            end
        end

        function set.computeAeroCoeffs(self, aeroModelFn)
            if Constants.VALIDATE_FLAG
                self.computeAeroCoeffs = Validator.validateType(aeroModelFn, "function_handle");
            else
                self.computeAeroCoeffs = aeroModelFn;
            end
        end
    end
end


% Note 1
%
% Parameter values are obtained by directly getting the value from the self.params vector using the
% respective parameter index. For example:
% 
% | m = self.params(self.mIdx);
% 
% This is extremely fast. A possible alternative is to abstract this direct indexing call behind a
% dependent "parameter" with a getter. For example, define:
%
% | properties (Dependent)
% |     m
% | end
% 
% | methods
% |     function m = get.m(self)
% |         m = self.params(self.mIdx);
% |     end
% | end
%
% Then, to get the parameter value:
%
% | m = self.m
%
% which calls the getter. This results is arguably more readable code, at the cost of a significant
% performance hit (~15% slower).
%
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
        nParams = 0;

        nEstimatedParams = 0;
        estimatedParamIdxs = [];

        dIdx = 0;
        SIdx = 0;
        nFinsIdx = 0;
        deltaFinsIdx = 0;

        mIdx = 0;
        IxxIdx = 0;

        CDIdx = 0;
        CDTable_Mach0Idx = 0;
        CDTable_CD0Idx = 0;
        CDTable_Len = 0;

        Cl0Idx = 0;
        Cl0Table_Mach0Idx = 0;
        Cl0Table_Cl0Idx = 0;
        Cl0Table_Len = 0;

        ClpIdx = 0;
        ClpTable_Mach0Idx = 0;
        ClpTable_ClpIdx = 0;
        ClpTable_Len = 0;

        CldeltaIdx = 0;
        CldeltaTable_Mach0Idx = 0;
        CldeltaTable_CldeltaIdx = 0;
        CldeltaTable_Len = 0;

        isAeroModelInit = false;

        computeAeroCoeffs
    end

    properties (Constant)
        nStates = 7;

        xIdx = 1;
        yIdx = 2;
        zIdx = 3;
        vxIdx = 4;
        vyIdx = 5;
        vzIdx = 6;
        pIdx = 7;

        DEFAULT_D = 0.02999232;
        DEFAULT_S = (pi / 4) * 0.02999232 ^ 2;
        DEFAULT_NFINS = 4;
        DEFAULT_DELTAFINS = 0;

        DEFAULT_M = 1.58885397;
        DEFAULT_IXX = 0.000192309;

        DEFAULT_CD = 0.472;
        DEFAULT_CD_TABLE_X = [0; 5];
        DEFAULT_CD_TABLE_Y = [0.472; 0.472];

        DEFAULT_Cl0 = 0;
        DEFAULT_Cl0_TABLE_X = [0; 5];
        DEFAULT_Cl0_TABLE_Y = [0; 0];

        DEFAULT_Clp = -4.5;
        DEFAULT_Clp_TABLE_X = [0; 5];
        DEFAULT_Clp_TABLE_Y = [-4.5; -4.5];

        DEFAULT_Cldelta = 0;
        DEFAULT_Cldelta_TABLE_X = [0; 5];
        DEFAULT_Cldelta_TABLE_Y = [0; 0];

        VALID_AERO_MODELS = ["constant", "table"];
    end
    

    methods
        % Constructor ==============================================================================

        function self = Projectile(aeroModel)
            self.stateDef = StateDef();

            self.paramDefs.d = ParamDef(self.DEFAULT_D);
            self.paramDefs.S = ParamDef(self.DEFAULT_S);
            self.paramDefs.nFins = ParamDef(self.DEFAULT_NFINS);
            self.paramDefs.deltaFins = ParamDef(self.DEFAULT_DELTAFINS);

            self.paramDefs.m = ParamDef(self.DEFAULT_M);
            self.paramDefs.Ixx = ParamDef(self.DEFAULT_IXX);
            
            if nargin == 0
                self.aeroModel = "constant";
            elseif nargin == 1
                self.aeroModel = aeroModel;
            else
                error("Too many input arguments.")
            end

            self.update();
        end

        
        % Update methods ===========================================================================

        function update(self)
            self.updateState();

            self.updateModels();
            self.updateParams();
            self.updateEstimatedParams();
        end


        function updateState(self)
            self.time = self.stateDef.time;
            self.state = self.stateDef.state;
            self.stateCovar = self.stateDef.covar;
        end

        
        function updateModels(self)
            self.updateAeroModel();
        end
        

        function updateAeroModel(self)
            switch self.aeroModel
                case "constant"
                    self.computeAeroCoeffs = @self.constantAeroModel;
                    
                    if ~self.isAeroModelInit
                        self.paramDefs.CD      = ParamDef(self.DEFAULT_CD);
                        self.paramDefs.Cl0     = ParamDef(self.DEFAULT_Cl0);
                        self.paramDefs.Clp     = ParamDef(self.DEFAULT_Clp);
                        self.paramDefs.Cldelta = ParamDef(self.DEFAULT_Cldelta);
                    end

                case "table"
                    self.computeAeroCoeffs = @self.tableAeroModel;
                    
                    if ~self.isAeroModelInit
                        self.paramDefs.CD      = ParamTableDef(self.DEFAULT_CD_TABLE_X,      self.DEFAULT_CD_TABLE_Y);
                        self.paramDefs.Cl0     = ParamTableDef(self.DEFAULT_Cl0_TABLE_X,     self.DEFAULT_Cl0_TABLE_Y);
                        self.paramDefs.Clp     = ParamTableDef(self.DEFAULT_Clp_TABLE_X,     self.DEFAULT_Clp_TABLE_Y);
                        self.paramDefs.Cldelta = ParamTableDef(self.DEFAULT_Cldelta_TABLE_X, self.DEFAULT_Cldelta_TABLE_Y);
                    end

                otherwise
                    error("Invalid aerodynamics model: %s.", self.aeroModel)
            end

            self.isAeroModelInit = true;
        end


        function updateParams(self)
            % See Note 1 regarding parameter indexing

            self.params = [];
            
            self.dIdx = 0;
            self.SIdx = 0;
            self.nFinsIdx = 0;
            self.deltaFinsIdx = 0;

            self.mIdx = 0;
            self.IxxIdx = 0;
    
            self.CDIdx = 0;
            self.CDTable_Mach0Idx = 0;
            self.CDTable_CD0Idx = 0;
            self.CDTable_Len = 0;

            self.Cl0Idx = 0;
            self.Cl0Table_Mach0Idx = 0;
            self.Cl0Table_Cl0Idx = 0;
            self.Cl0Table_Len = 0;
    
            self.ClpIdx = 0;
            self.ClpTable_Mach0Idx = 0;
            self.ClpTable_ClpIdx = 0;
            self.ClpTable_Len = 0;
    
            self.CldeltaIdx = 0;
            self.CldeltaTable_Mach0Idx = 0;
            self.CldeltaTable_CldeltaIdx = 0;
            self.CldeltaTable_Len = 0;
            
            % d
            self.dIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.d.value];
            
            % S
            self.SIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.S.value];

            % nFins
            self.nFinsIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.nFins.value];

            % deltaFins
            self.deltaFinsIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.deltaFins.value];

            % m
            self.mIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.m.value];

            % I
            self.IxxIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.Ixx.value];
            
            switch self.aeroModel
                case "constant"
                    % CD
                    self.CDIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CD.value];

                    % Cl0
                    self.Cl0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cl0.value];

                    % Clp
                    self.ClpIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Clp.value];

                    % Cldelta
                    self.CldeltaIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cldelta.value];

                case "table"
                    % CD
                    self.CDTable_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CD.xValues];

                    self.CDTable_CD0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CD.yValues];
                    
                    self.CDTable_Len = self.paramDefs.CD.nValues;

                    % Cl0
                    self.Cl0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cl0.xValues];

                    self.Cl0Table_Cl0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cl0.yValues];
                    
                    self.Cl0Table_Len = self.paramDefs.Cl0.nValues;

                    % Clp
                    self.ClpTable_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Clp.xValues];

                    self.ClpTable_ClpIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Clp.yValues];
                    
                    self.ClpTable_Len = self.paramDefs.Clp.nValues;

                    % Cldelta
                    self.CldeltaTable_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cldelta.xValues];

                    self.CldeltaTable_CldeltaIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cldelta.yValues];
                    
                    self.CldeltaTable_Len = self.paramDefs.Cldelta.nValues;
            end
        end


        function updateEstimatedParams(self)
            self.estimatedParams = [];
            self.estimatedParamCovar = [];
            self.estimatedParamIdxs = [];

            % d
            if self.paramDefs.d.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.d.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.d.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.dIdx];
            end

            % S
            if self.paramDefs.S.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.S.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.S.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.SIdx];
            end

            % nFins
            if self.paramDefs.nFins.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.nFins.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.nFins.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.nFinsIdx];
            end

            % deltaFins
            if self.paramDefs.deltaFins.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.deltaFins.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.deltaFins.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.deltaFinsIdx];
            end

            % m
            if self.paramDefs.m.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.m.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.m.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.mIdx];
            end

            % I
            if self.paramDefs.Ixx.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.Ixx.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Ixx.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.IxxIdx];
            end
            
            switch self.aeroModel
                case "constant"
                    % CD
                    if self.paramDefs.CD.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CD.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CD.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CDIdx];
                    end

                    % Cl0
                    if self.paramDefs.Cl0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.Cl0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cl0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.Cl0Idx];
                    end

                    % Clp
                    if self.paramDefs.Clp.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.Clp.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Clp.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.ClpIdx];
                    end

                    % Cldelta
                    if self.paramDefs.Cldelta.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.Cldelta.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cldelta.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CldeltaIdx];
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

                    % Cl0
                    for i = 1:length(self.paramDefs.Cl0.yValues)
                        if self.paramDefs.Cl0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.Cl0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cl0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.Cl0Table_Cl0Idx + (i - 1)];
                        end
                    end

                    % Clp
                    for i = 1:length(self.paramDefs.Clp.yValues)
                        if self.paramDefs.Clp.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.Clp.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Clp.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.ClpTable_ClpIdx + (i - 1)];
                        end
                    end

                    % Cldelta
                    for i = 1:length(self.paramDefs.Cldelta.yValues)
                        if self.paramDefs.Cldelta.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.Cldelta.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cldelta.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CldeltaTable_CldeltaIdx + (i - 1)];
                        end
                    end
            end
        end


        function readPropsFromFile(self, filePath)
            try
                props = readmatrix(filePath);
                props = props(:, 2);  % Remove header column
            catch
                error("Cannot read properties CSV file. File either does not exist or is not formatted properly.")
            end

            self.paramDefs.d         = ParamDef(props(1));
            self.paramDefs.S         = ParamDef((pi / 4) * props(1) ^ 2);
            self.paramDefs.nFins     = ParamDef(props(2));
            self.paramDefs.deltaFins = ParamDef(props(3));
            self.paramDefs.m         = ParamDef(props(4));
            self.paramDefs.Ixx       = ParamDef(props(5));
        end


        function readAeroModelTablesFromFile(self, filePath)
            try
                aeroTable = readmatrix(filePath);
            catch
                error("Cannot read aero table CSV file. File either does not exist or is not formatted properly.")
            end

            machValues    =  aeroTable(:, 1);
            CDValues      = -aeroTable(:, 2);
            Cl0Values     =  aeroTable(:, 10);
            ClpValues     =  aeroTable(:, 11);
            CldeltaValues =  aeroTable(:, 12);

            self.paramDefs.CD      = ParamTableDef(machValues, CDValues);
            self.paramDefs.Cl0     = ParamTableDef(machValues, Cl0Values);
            self.paramDefs.Clp     = ParamTableDef(machValues, ClpValues);
            self.paramDefs.Cldelta = ParamTableDef(machValues, CldeltaValues);
        end

        
        % Model methods ============================================================================

        function [CD, Cl0, Clp, Cldelta] = constantAeroModel(self, ~)
            CD = self.params(self.CDIdx);
            Cl0 = self.params(self.Cl0Idx);
            Clp = self.params(self.ClpIdx);
            Cldelta = self.params(self.CldeltaIdx);
        end


        function [CD, Cl0, Clp, Cldelta] = tableAeroModel(self, mach)
            dmach = self.params(self.CDTable_Mach0Idx + 1) - self.params(self.CDTable_Mach0Idx);  % TODO: Non-uniform table

            CD = 0;
            for i = 0:(self.CDTable_Len - 1)
                mach_i = self.params(self.CDTable_Mach0Idx + i);
                CD_i = self.params(self.CDTable_CD0Idx + i);

                CD = CD + CD_i * self.linearKernel((mach - mach_i) / dmach);
            end

            Cl0 = 0;
            for i = 0:(self.Cl0Table_Len - 1)
                mach_i = self.params(self.Cl0Table_Mach0Idx + i);
                Cl0_i = self.params(self.Cl0Table_Cl0Idx + i);

                Cl0 = Cl0 + Cl0_i * self.linearKernel((mach - mach_i) / dmach);
            end

            Clp = 0;
            for i = 0:(self.ClpTable_Len - 1)
                mach_i = self.params(self.ClpTable_Mach0Idx + i);
                Clp_i = self.params(self.ClpTable_ClpIdx + i);

                Clp = Clp + Clp_i * self.linearKernel((mach - mach_i) / dmach);
            end

            Cldelta = 0;
            for i = 0:(self.CldeltaTable_Len - 1)
                mach_i = self.params(self.CldeltaTable_Mach0Idx + i);
                Cldelta_i = self.params(self.CldeltaTable_CldeltaIdx + i);

                Cldelta = Cldelta + Cldelta_i * self.linearKernel((mach - mach_i) / dmach);
            end
        end
        

        function k = linearKernel(~, x)
            if (-1 < x) && (x < 0)
                k = 1 + x;
            elseif (0 <= x) && (x < 1)
                k = 1 - x;
            else
                k = 0;
            end
        end

        
        % Setters ==================================================================================

        function set.stateDef(self, stateDef)
            if Settings.VALIDATE_FLAG
                self.stateDef = Validator.validateType(stateDef, "StateDef");
            else
                self.stateDef = stateDef;
            end
        end

        function set.paramDefs(self, paramDefs)
            if Settings.VALIDATE_FLAG
                self.paramDefs = Validator.validateFieldTypes(paramDefs, ["ParamDef", "ParamTableDef"]);
            else
                self.paramDefs = paramDefs;
            end
        end

        function set.state(self, state)
            if Settings.VALIDATE_FLAG
                state = Validator.validateType(state, "double");
                self.state = Validator.validateSize(state, [self.nStates, 1]);
            else
                self.state = state;
            end
        end

        function set.stateCovar(self, stateCovar)
            if Settings.VALIDATE_FLAG
                stateCovar = Validator.validateType(stateCovar, "double");
                self.stateCovar = Validator.validateSize(stateCovar, [self.nStates, self.nStates]);
            else
                self.stateCovar = stateCovar;
            end
        end

        function set.params(self, params)
            if Settings.VALIDATE_FLAG
                self.params = Validator.validateType(params, "double");
            else
                self.params = params;
            end

            self.nParams = length(params);
        end

        function set.estimatedParams(self, estimatedParams)
            if Settings.VALIDATE_FLAG
                self.estimatedParams = Validator.validateType(estimatedParams, "double");
            else
                self.estimatedParams = estimatedParams;
            end

            self.nEstimatedParams = length(estimatedParams);
        end

        function set.estimatedParamCovar(self, estimatedParamCovar)
            if Settings.VALIDATE_FLAG
                estimatedParamCovar = Validator.validateType(estimatedParamCovar, "double");
                self.estimatedParamCovar = Validator.validateSize(estimatedParamCovar, [self.nEstimatedParams, self.nEstimatedParams]);
            else
                self.estimatedParamCovar = estimatedParamCovar;
            end
        end

        function set.estimatedParamIdxs(self, estimatedParamIdxs)
            if Settings.VALIDATE_FLAG
                self.estimatedParamIdxs = Validator.validateType(estimatedParamIdxs, "double");
            else
                self.estimatedParamIdxs = estimatedParamIdxs;
            end
        end

        function set.aeroModel(self, aeroModel)
            if Settings.VALIDATE_FLAG
                self.aeroModel = Validator.validateString(aeroModel, self.VALID_AERO_MODELS);
            else
                self.aeroModel = aeroModel;
            end

            self.isAeroModelInit = false;
            self.updateAeroModel();
        end

        function set.computeAeroCoeffs(self, aeroModelFn)
            if Settings.VALIDATE_FLAG
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
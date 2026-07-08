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

        consideredParams = [];
        consideredParamCovar = [];

        aeroModel
    end

    properties (SetAccess = private)
        nParams = 0;

        nEstimatedParams = 0;
        estimatedParamIdxs = [];

        nConsideredParams = 0;
        consideredParamIdxs = [];
        
        dIdx = 0;
        SIdx = 0;
        nFinsIdx = 0;
        deltaFinsIdx = 0;

        mIdx = 0;
        IxxIdx = 0;
        IyyIdx = 0;
        IzzIdx = 0;
        IxyIdx = 0;
        IxzIdx = 0;
        IyzIdx = 0;

        CX0Idx = 0;
        CX0Table_Mach0Idx = 0;
        CX0Table_CX0Idx = 0;
        CX0Table_Len = 0;

        CX2Idx = 0;
        CX2Table_Mach0Idx = 0;
        CX2Table_CX2Idx = 0;
        CX2Table_Len = 0;

        CY0Idx = 0;
        CY0Table_Mach0Idx = 0;
        CY0Table_CY0Idx = 0;
        CY0Table_Len = 0;

        CZ0Idx = 0;
        CZ0Table_Mach0Idx = 0;
        CZ0Table_CZ0Idx = 0;
        CZ0Table_Len = 0;

        CNalpha0Idx = 0;
        CNalpha0Table_Mach0Idx = 0;
        CNalpha0Table_CNalpha0Idx = 0;
        CNalpha0Table_Len = 0;

        CNalpha2Idx = 0;
        CNalpha2Table_Mach0Idx = 0;
        CNalpha2Table_CNalpha2Idx = 0;
        CNalpha2Table_Len = 0;

        CNpalpha0Idx = 0;
        CNpalpha0Table_Mach0Idx = 0;
        CNpalpha0Table_CNpalpha0Idx = 0;
        CNpalpha0Table_Len = 0;

        CNpalpha2Idx = 0;
        CNpalpha2Table_Mach0Idx = 0;
        CNpalpha2Table_CNpalpha2Idx = 0;
        CNpalpha2Table_Len = 0;

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

        Cm0Idx = 0;
        Cm0Table_Mach0Idx = 0;
        Cm0Table_Cm0Idx = 0;
        Cm0Table_Len = 0;

        Cn0Idx = 0;
        Cn0Table_Mach0Idx = 0;
        Cn0Table_Cn0Idx = 0;
        Cn0Table_Len = 0;

        CMalpha0Idx = 0;
        CMalpha0Table_Mach0Idx = 0;
        CMalpha0Table_CMalpha0Idx = 0;
        CMalpha0Table_Len = 0;

        CMalpha2Idx = 0;
        CMalpha2Table_Mach0Idx = 0;
        CMalpha2Table_CMalpha2Idx = 0;
        CMalpha2Table_Len = 0;

        CMpalpha0Idx = 0;
        CMpalpha0Table_Mach0Idx = 0;
        CMpalpha0Table_CMpalpha0Idx = 0;
        CMpalpha0Table_Len = 0;

        CMpalpha2Idx = 0;
        CMpalpha2Table_Mach0Idx = 0;
        CMpalpha2Table_CMpalpha2Idx = 0;
        CMpalpha2Table_Len = 0;

        CMqIdx = 0;
        CMqTable_Mach0Idx = 0;
        CMqTable_CMqIdx = 0;
        CMqTable_Len = 0;

        isAeroModelInit = false;

        computeAeroCoeffs
    end

    properties (Constant)
        nStates = 12;

        xIdx = 1;
        yIdx = 2;
        zIdx = 3;
        phiIdx = 4;
        thetaIdx = 5;
        psiIdx = 6;
        vxIdx = 7;
        vyIdx = 8;
        vzIdx = 9;
        pIdx = 10;
        qIdx = 11;
        rIdx = 12;

        DEFAULT_D = 0.02999232;
        DEFAULT_S = (pi / 4) * 0.02999232 ^ 2;
        DEFAULT_NFINS = 4;
        DEFAULT_DELTAFINS = 0;

        DEFAULT_M = 1.58885397;
        DEFAULT_IXX = 0.000192309;
        DEFAULT_IYY = 0.00986927;
        DEFAULT_IZZ = 0.00986927;
        DEFAULT_IXY = 0;
        DEFAULT_IXZ = 0;
        DEFAULT_IYZ = 0;

        DEFAULT_CX0 = -0.472;
        DEFAULT_CX0_TABLE_X = [0; 5];
        DEFAULT_CX0_TABLE_Y = [-0.472; -0.472];

        DEFAULT_CX2 = -3.32;
        DEFAULT_CX2_TABLE_X = [0; 5];
        DEFAULT_CX2_TABLE_Y = [-3.32; -3.32];

        DEFAULT_CY0 = 0;
        DEFAULT_CY0_TABLE_X = [0; 5];
        DEFAULT_CY0_TABLE_Y = [0; 0];

        DEFAULT_CZ0 = 0;
        DEFAULT_CZ0_TABLE_X = [0; 5];
        DEFAULT_CZ0_TABLE_Y = [0; 0];

        DEFAULT_CNalpha0 = 13.71;
        DEFAULT_CNalpha0_TABLE_X = [0; 5];
        DEFAULT_CNalpha0_TABLE_Y = [13.71; 13.71];

        DEFAULT_CNalpha2 = 0;
        DEFAULT_CNalpha2_TABLE_X = [0; 5];
        DEFAULT_CNalpha2_TABLE_Y = [0; 0];

        DEFAULT_CNpalpha0 = 0;
        DEFAULT_CNpalpha0_TABLE_X = [0; 5];
        DEFAULT_CNpalpha0_TABLE_Y = [0; 0];

        DEFAULT_CNpalpha2 = 0;
        DEFAULT_CNpalpha2_TABLE_X = [0; 5];
        DEFAULT_CNpalpha2_TABLE_Y = [0; 0];

        DEFAULT_Cl0 = 0;
        DEFAULT_Cl0_TABLE_X = [0; 5];
        DEFAULT_Cl0_TABLE_Y = [0; 0];

        DEFAULT_Clp = -4.5;
        DEFAULT_Clp_TABLE_X = [0; 5];
        DEFAULT_Clp_TABLE_Y = [-4.5; -4.5];

        DEFAULT_Cldelta = 0;
        DEFAULT_Cldelta_TABLE_X = [0; 5];
        DEFAULT_Cldelta_TABLE_Y = [0; 0];

        DEFAULT_Cm0 = 0;
        DEFAULT_Cm0_TABLE_X = [0; 5];
        DEFAULT_Cm0_TABLE_Y = [0; 0];

        DEFAULT_Cn0 = 0;
        DEFAULT_Cn0_TABLE_X = [0; 5];
        DEFAULT_Cn0_TABLE_Y = [0; 0];

        DEFAULT_CMalpha0 = -22.01381098;
        DEFAULT_CMalpha0_TABLE_X = [0; 5];
        DEFAULT_CMalpha0_TABLE_Y = [-22.01381098; -22.01381098];

        DEFAULT_CMalpha2 = 0;
        DEFAULT_CMalpha2_TABLE_X = [0; 5];
        DEFAULT_CMalpha2_TABLE_Y = [0; 0];

        DEFAULT_CMpalpha0 = 0;
        DEFAULT_CMpalpha0_TABLE_X = [0; 5];
        DEFAULT_CMpalpha0_TABLE_Y = [0; 0];

        DEFAULT_CMpalpha2 = 0;
        DEFAULT_CMpalpha2_TABLE_X = [0; 5];
        DEFAULT_CMpalpha2_TABLE_Y = [0; 0];

        DEFAULT_CMq = -207.1;
        DEFAULT_CMq_TABLE_X = [0; 5];
        DEFAULT_CMq_TABLE_Y = [-207.1; -207.1];

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
            self.paramDefs.Iyy = ParamDef(self.DEFAULT_IYY);
            self.paramDefs.Izz = ParamDef(self.DEFAULT_IZZ);
            self.paramDefs.Ixy = ParamDef(self.DEFAULT_IXY);
            self.paramDefs.Ixz = ParamDef(self.DEFAULT_IXZ);
            self.paramDefs.Iyz = ParamDef(self.DEFAULT_IYZ);
            
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
            self.updateConsideredParams();
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
                        self.paramDefs.CX0       = ParamDef(self.DEFAULT_CX0);
                        self.paramDefs.CX2       = ParamDef(self.DEFAULT_CX2);
                        self.paramDefs.CY0       = ParamDef(self.DEFAULT_CY0);
                        self.paramDefs.CZ0       = ParamDef(self.DEFAULT_CZ0);
                        self.paramDefs.CNalpha0  = ParamDef(self.DEFAULT_CNalpha0);
                        self.paramDefs.CNalpha2  = ParamDef(self.DEFAULT_CNalpha2);
                        self.paramDefs.CNpalpha0 = ParamDef(self.DEFAULT_CNpalpha0);
                        self.paramDefs.CNpalpha2 = ParamDef(self.DEFAULT_CNpalpha2);
                        self.paramDefs.Cl0       = ParamDef(self.DEFAULT_Cl0);
                        self.paramDefs.Clp       = ParamDef(self.DEFAULT_Clp);
                        self.paramDefs.Cldelta   = ParamDef(self.DEFAULT_Cldelta);
                        self.paramDefs.Cm0       = ParamDef(self.DEFAULT_Cm0);
                        self.paramDefs.Cn0       = ParamDef(self.DEFAULT_Cn0);
                        self.paramDefs.CMalpha0  = ParamDef(self.DEFAULT_CMalpha0);
                        self.paramDefs.CMalpha2  = ParamDef(self.DEFAULT_CMalpha2);
                        self.paramDefs.CMpalpha0 = ParamDef(self.DEFAULT_CMpalpha0);
                        self.paramDefs.CMpalpha2 = ParamDef(self.DEFAULT_CMpalpha2);
                        self.paramDefs.CMq       = ParamDef(self.DEFAULT_CMq);
                    end

                case "table"
                    self.computeAeroCoeffs = @self.tableAeroModel;
                    
                    if ~self.isAeroModelInit
                        self.paramDefs.CX0       = ParamTableDef(self.DEFAULT_CX0_TABLE_X,       self.DEFAULT_CX0_TABLE_Y);
                        self.paramDefs.CX2       = ParamTableDef(self.DEFAULT_CX2_TABLE_X,       self.DEFAULT_CX2_TABLE_Y);
                        self.paramDefs.CY0       = ParamTableDef(self.DEFAULT_CY0_TABLE_X,       self.DEFAULT_CY0_TABLE_Y);
                        self.paramDefs.CZ0       = ParamTableDef(self.DEFAULT_CZ0_TABLE_X,       self.DEFAULT_CZ0_TABLE_Y);
                        self.paramDefs.CNalpha0  = ParamTableDef(self.DEFAULT_CNalpha0_TABLE_X,  self.DEFAULT_CNalpha0_TABLE_Y);
                        self.paramDefs.CNalpha2  = ParamTableDef(self.DEFAULT_CNalpha2_TABLE_X,  self.DEFAULT_CNalpha2_TABLE_Y);
                        self.paramDefs.CNpalpha0 = ParamTableDef(self.DEFAULT_CNpalpha0_TABLE_X, self.DEFAULT_CNpalpha0_TABLE_Y);
                        self.paramDefs.CNpalpha2 = ParamTableDef(self.DEFAULT_CNpalpha2_TABLE_X, self.DEFAULT_CNpalpha2_TABLE_Y);
                        self.paramDefs.Cl0       = ParamTableDef(self.DEFAULT_Cl0_TABLE_X,       self.DEFAULT_Cl0_TABLE_Y);
                        self.paramDefs.Clp       = ParamTableDef(self.DEFAULT_Clp_TABLE_X,       self.DEFAULT_Clp_TABLE_Y);
                        self.paramDefs.Cldelta   = ParamTableDef(self.DEFAULT_Cldelta_TABLE_X,   self.DEFAULT_Cldelta_TABLE_Y);
                        self.paramDefs.Cm0       = ParamTableDef(self.DEFAULT_Cm0_TABLE_X,       self.DEFAULT_Cm0_TABLE_Y);
                        self.paramDefs.Cn0       = ParamTableDef(self.DEFAULT_Cn0_TABLE_X,       self.DEFAULT_Cn0_TABLE_Y);
                        self.paramDefs.CMalpha0  = ParamTableDef(self.DEFAULT_CMalpha0_TABLE_X,  self.DEFAULT_CMalpha0_TABLE_Y);
                        self.paramDefs.CMalpha2  = ParamTableDef(self.DEFAULT_CMalpha2_TABLE_X,  self.DEFAULT_CMalpha2_TABLE_Y);
                        self.paramDefs.CMpalpha0 = ParamTableDef(self.DEFAULT_CMpalpha0_TABLE_X, self.DEFAULT_CMpalpha0_TABLE_Y);
                        self.paramDefs.CMpalpha2 = ParamTableDef(self.DEFAULT_CMpalpha2_TABLE_X, self.DEFAULT_CMpalpha2_TABLE_Y);
                        self.paramDefs.CMq       = ParamTableDef(self.DEFAULT_CMq_TABLE_X,       self.DEFAULT_CMq_TABLE_Y);
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
            self.IyyIdx = 0;
            self.IzzIdx = 0;
            self.IxyIdx = 0;
            self.IxzIdx = 0;
            self.IyzIdx = 0;
    
            self.CX0Idx = 0;
            self.CX0Table_Mach0Idx = 0;
            self.CX0Table_CX0Idx = 0;
            self.CX0Table_Len = 0;

            self.CX2Idx = 0;
            self.CX2Table_Mach0Idx = 0;
            self.CX2Table_CX2Idx = 0;
            self.CX2Table_Len = 0;
    
            self.CY0Idx = 0;
            self.CY0Table_Mach0Idx = 0;
            self.CY0Table_CY0Idx = 0;
            self.CY0Table_Len = 0;
    
            self.CZ0Idx = 0;
            self.CZ0Table_Mach0Idx = 0;
            self.CZ0Table_CZ0Idx = 0;
            self.CZ0Table_Len = 0;
    
            self.CNalpha0Idx = 0;
            self.CNalpha0Table_Mach0Idx = 0;
            self.CNalpha0Table_CNalpha0Idx = 0;
            self.CNalpha0Table_Len = 0;
    
            self.CNalpha2Idx = 0;
            self.CNalpha2Table_Mach0Idx = 0;
            self.CNalpha2Table_CNalpha2Idx = 0;
            self.CNalpha2Table_Len = 0;
    
            self.CNpalpha0Idx = 0;
            self.CNpalpha0Table_Mach0Idx = 0;
            self.CNpalpha0Table_CNpalpha0Idx = 0;
            self.CNpalpha0Table_Len = 0;
    
            self.CNpalpha2Idx = 0;
            self.CNpalpha2Table_Mach0Idx = 0;
            self.CNpalpha2Table_CNpalpha2Idx = 0;
            self.CNpalpha2Table_Len = 0;
    
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
    
            self.Cm0Idx = 0;
            self.Cm0Table_Mach0Idx = 0;
            self.Cm0Table_Cm0Idx = 0;
            self.Cm0Table_Len = 0;
    
            self.Cn0Idx = 0;
            self.Cn0Table_Mach0Idx = 0;
            self.Cn0Table_Cn0Idx = 0;
            self.Cn0Table_Len = 0;
    
            self.CMalpha0Idx = 0;
            self.CMalpha0Table_Mach0Idx = 0;
            self.CMalpha0Table_CMalpha0Idx = 0;
            self.CMalpha0Table_Len = 0;
    
            self.CMalpha2Idx = 0;
            self.CMalpha2Table_Mach0Idx = 0;
            self.CMalpha2Table_CMalpha2Idx = 0;
            self.CMalpha2Table_Len = 0;
    
            self.CMpalpha0Idx = 0;
            self.CMpalpha0Table_Mach0Idx = 0;
            self.CMpalpha0Table_CMpalpha0Idx = 0;
            self.CMpalpha0Table_Len = 0;
    
            self.CMpalpha2Idx = 0;
            self.CMpalpha2Table_Mach0Idx = 0;
            self.CMpalpha2Table_CMpalpha2Idx = 0;
            self.CMpalpha2Table_Len = 0;
    
            self.CMqIdx = 0;
            self.CMqTable_Mach0Idx = 0;
            self.CMqTable_CMqIdx = 0;
            self.CMqTable_Len = 0;

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

            self.IyyIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.Iyy.value];

            self.IzzIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.Izz.value];

            self.IxyIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.Ixy.value];

            self.IxzIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.Ixz.value];

            self.IyzIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.Iyz.value];
            
            switch self.aeroModel
                case "constant"
                    % CX0
                    self.CX0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CX0.value];

                    % CX2
                    self.CX2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CX2.value];
                    
                    % CY0
                    self.CY0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CY0.value];

                    % CZ0
                    self.CZ0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CZ0.value];

                    % CNalpha0
                    self.CNalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNalpha0.value];

                    % CNalpha2
                    self.CNalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNalpha2.value];

                    % CNpalpha0
                    self.CNpalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNpalpha0.value];

                    % CNpalpha2
                    self.CNpalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNpalpha2.value];

                    % Cl0
                    self.Cl0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cl0.value];

                    % Clp
                    self.ClpIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Clp.value];

                    % Cldelta
                    self.CldeltaIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cldelta.value];

                    % Cm0
                    self.Cm0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cm0.value];

                    % Cn0
                    self.Cn0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cn0.value];

                    % CMalpha0
                    self.CMalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMalpha0.value];

                    % CMalpha2
                    self.CMalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMalpha2.value];

                    % CMpalpha0
                    self.CMpalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMpalpha0.value];

                    % CMpalpha2
                    self.CMpalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMpalpha2.value];

                    % CMq
                    self.CMqIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMq.value];

                case "table"
                    % CX0
                    self.CX0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CX0.xValues];

                    self.CX0Table_CX0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CX0.yValues];
                    
                    self.CX0Table_Len = self.paramDefs.CX0.nValues;

                    % CX2
                    self.CX2Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CX2.xValues];

                    self.CX2Table_CX2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CX2.yValues];
                    
                    self.CX2Table_Len = self.paramDefs.CX2.nValues;

                    % CY0
                    self.CY0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CY0.xValues];

                    self.CY0Table_CY0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CY0.yValues];
                    
                    self.CY0Table_Len = self.paramDefs.CY0.nValues;

                    % CZ0
                    self.CZ0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CZ0.xValues];

                    self.CZ0Table_CZ0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CZ0.yValues];
                    
                    self.CZ0Table_Len = self.paramDefs.CZ0.nValues;

                    % CNalpha0
                    self.CNalpha0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNalpha0.xValues];

                    self.CNalpha0Table_CNalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNalpha0.yValues];
                    
                    self.CNalpha0Table_Len = self.paramDefs.CNalpha0.nValues;

                    % CNalpha2
                    self.CNalpha2Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNalpha2.xValues];

                    self.CNalpha2Table_CNalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNalpha2.yValues];
                    
                    self.CNalpha2Table_Len = self.paramDefs.CNalpha2.nValues;

                    % CNpalpha0
                    self.CNpalpha0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNpalpha0.xValues];

                    self.CNpalpha0Table_CNpalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNpalpha0.yValues];
                    
                    self.CNpalpha0Table_Len = self.paramDefs.CNpalpha0.nValues;

                    % CNpalpha2
                    self.CNpalpha2Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNpalpha2.xValues];

                    self.CNpalpha2Table_CNpalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CNpalpha2.yValues];
                    
                    self.CNpalpha2Table_Len = self.paramDefs.CNpalpha2.nValues;

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

                    % Cm0
                    self.Cm0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cm0.xValues];

                    self.Cm0Table_Cm0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cm0.yValues];
                    
                    self.Cm0Table_Len = self.paramDefs.Cm0.nValues;

                    % Cn0
                    self.Cn0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cn0.xValues];

                    self.Cn0Table_Cn0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.Cn0.yValues];
                    
                    self.Cn0Table_Len = self.paramDefs.Cn0.nValues;

                    % CMalpha0
                    self.CMalpha0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMalpha0.xValues];

                    self.CMalpha0Table_CMalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMalpha0.yValues];
                    
                    self.CMalpha0Table_Len = self.paramDefs.CMalpha0.nValues;

                    % CMalpha2
                    self.CMalpha2Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMalpha2.xValues];

                    self.CMalpha2Table_CMalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMalpha2.yValues];
                    
                    self.CMalpha2Table_Len = self.paramDefs.CMalpha2.nValues;

                    % CMpalpha0
                    self.CMpalpha0Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMpalpha0.xValues];

                    self.CMpalpha0Table_CMpalpha0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMpalpha0.yValues];
                    
                    self.CMpalpha0Table_Len = self.paramDefs.CMpalpha0.nValues;

                    % CMpalpha2
                    self.CMpalpha2Table_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMpalpha2.xValues];

                    self.CMpalpha2Table_CMpalpha2Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMpalpha2.yValues];
                    
                    self.CMpalpha2Table_Len = self.paramDefs.CMpalpha2.nValues;

                    % CMq
                    self.CMqTable_Mach0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMq.xValues];

                    self.CMqTable_CMqIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.CMq.yValues];
                    
                    self.CMqTable_Len = self.paramDefs.CMq.nValues;
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

            if self.paramDefs.Iyy.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.Iyy.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Iyy.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.IyyIdx];
            end

            if self.paramDefs.Izz.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.Izz.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Izz.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.IzzIdx];
            end

            if self.paramDefs.Ixy.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.Ixy.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Ixy.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.IxyIdx];
            end

            if self.paramDefs.Ixz.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.Ixz.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Ixz.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.IxzIdx];
            end

            if self.paramDefs.Iyz.isEstimated
                self.estimatedParams = [self.estimatedParams; self.paramDefs.Iyz.value];
                self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Iyz.covar);
                self.estimatedParamIdxs = [self.estimatedParamIdxs; self.IyzIdx];
            end
            
            switch self.aeroModel
                case "constant"
                    % CX0
                    if self.paramDefs.CX0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CX0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CX0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CX0Idx];
                    end

                    % CX2
                    if self.paramDefs.CX2.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CX2.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CX2.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CX2Idx];
                    end

                    % CY0
                    if self.paramDefs.CY0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CY0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CY0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CY0Idx];
                    end

                    % CZ0
                    if self.paramDefs.CZ0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CZ0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CZ0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CZ0Idx];
                    end

                    % CNalpha0
                    if self.paramDefs.CNalpha0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CNalpha0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNalpha0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNalpha0Idx];
                    end

                    % CNalpha2
                    if self.paramDefs.CNalpha2.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CNalpha2.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNalpha2.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNalpha2Idx];
                    end

                    % CNpalpha0
                    if self.paramDefs.CNpalpha0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CNpalpha0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNpalpha0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNpalpha0Idx];
                    end

                    % CNpalpha2
                    if self.paramDefs.CNpalpha2.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CNpalpha2.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNpalpha2.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNpalpha2Idx];
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

                    % Cm0
                    if self.paramDefs.Cm0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.Cm0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cm0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.Cm0Idx];
                    end

                    % Cn0
                    if self.paramDefs.Cn0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.Cn0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cn0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.Cn0Idx];
                    end

                    % CMalpha0
                    if self.paramDefs.CMalpha0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CMalpha0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMalpha0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMalpha0Idx];
                    end

                    % CMalpha2
                    if self.paramDefs.CMalpha2.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CMalpha2.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMalpha2.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMalpha2Idx];
                    end

                    % CMpalpha0
                    if self.paramDefs.CMpalpha0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CMpalpha0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMpalpha0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMpalpha0Idx];
                    end

                    % CMpalpha2
                    if self.paramDefs.CMpalpha2.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CMpalpha2.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMpalpha2.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMpalpha2Idx];
                    end

                    % CMq
                    if self.paramDefs.CMq.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.CMq.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMq.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMqIdx];
                    end

                case "table"
                    % CX0
                    for i = 1:length(self.paramDefs.CX0.yValues)
                        if self.paramDefs.CX0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CX0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CX0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CX0Table_CX0Idx + (i - 1)];
                        end
                    end

                    % CX2
                    for i = 1:length(self.paramDefs.CX2.yValues)
                        if self.paramDefs.CX2.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CX2.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CX2.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CX2Table_CX2Idx + (i - 1)];
                        end
                    end

                    % CY0
                    for i = 1:length(self.paramDefs.CY0.yValues)
                        if self.paramDefs.CY0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CY0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CY0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CY0Table_CY0Idx + (i - 1)];
                        end
                    end

                    % CZ0
                    for i = 1:length(self.paramDefs.CZ0.yValues)
                        if self.paramDefs.CZ0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CZ0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CZ0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CZ0Table_CZ0Idx + (i - 1)];
                        end
                    end

                    % CNalpha0
                    for i = 1:length(self.paramDefs.CNalpha0.yValues)
                        if self.paramDefs.CNalpha0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CNalpha0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNalpha0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNalpha0Table_CNalpha0Idx + (i - 1)];
                        end
                    end

                    % CNalpha2
                    for i = 1:length(self.paramDefs.CNalpha2.yValues)
                        if self.paramDefs.CNalpha2.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CNalpha2.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNalpha2.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNalpha2Table_CNalpha2Idx + (i - 1)];
                        end
                    end

                    % CNpalpha0
                    for i = 1:length(self.paramDefs.CNpalpha0.yValues)
                        if self.paramDefs.CNpalpha0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CNpalpha0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNpalpha0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNpalpha0Table_CNpalpha0Idx + (i - 1)];
                        end
                    end

                    % CNpalpha2
                    for i = 1:length(self.paramDefs.CNpalpha2.yValues)
                        if self.paramDefs.CNpalpha2.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CNpalpha2.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CNpalpha2.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CNpalpha2Table_CNpalpha2Idx + (i - 1)];
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

                    % Cm0
                    for i = 1:length(self.paramDefs.Cm0.yValues)
                        if self.paramDefs.Cm0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.Cm0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cm0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.Cm0Table_Cm0Idx + (i - 1)];
                        end
                    end

                    % Cn0
                    for i = 1:length(self.paramDefs.Cn0.yValues)
                        if self.paramDefs.Cn0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.Cn0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.Cn0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.Cn0Table_Cn0Idx + (i - 1)];
                        end
                    end

                    % CMalpha0
                    for i = 1:length(self.paramDefs.CMalpha0.yValues)
                        if self.paramDefs.CMalpha0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CMalpha0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMalpha0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMalpha0Table_CMalpha0Idx + (i - 1)];
                        end
                    end

                    % CMalpha2
                    for i = 1:length(self.paramDefs.CMalpha2.yValues)
                        if self.paramDefs.CMalpha2.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CMalpha2.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMalpha2.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMalpha2Table_CMalpha2Idx + (i - 1)];
                        end
                    end

                    % CMpalpha0
                    for i = 1:length(self.paramDefs.CMpalpha0.yValues)
                        if self.paramDefs.CMpalpha0.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CMpalpha0.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMpalpha0.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMpalpha0Table_CMpalpha0Idx + (i - 1)];
                        end
                    end

                    % CMpalpha2
                    for i = 1:length(self.paramDefs.CMpalpha2.yValues)
                        if self.paramDefs.CMpalpha2.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CMpalpha2.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMpalpha2.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMpalpha2Table_CMpalpha2Idx + (i - 1)];
                        end
                    end

                    % CMq
                    for i = 1:length(self.paramDefs.CMq.yValues)
                        if self.paramDefs.CMq.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.CMq.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.CMq.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.CMqTable_CMqIdx + (i - 1)];
                        end
                    end
            end
        end


        function updateConsideredParams(self)
            self.consideredParams = [];
            self.consideredParamCovar = [];
            self.consideredParamIdxs = [];
            
            % d
            if self.paramDefs.d.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.d.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.d.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.dIdx];
            end

            % S
            if self.paramDefs.S.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.S.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.S.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.SIdx];
            end

            % nFins
            if self.paramDefs.nFins.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.nFins.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.nFins.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.nFinsIdx];
            end

            % deltaFins
            if self.paramDefs.deltaFins.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.deltaFins.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.deltaFins.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.deltaFinsIdx];
            end

            % m
            if self.paramDefs.m.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.m.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.m.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.mIdx];
            end

            % I
            if self.paramDefs.Ixx.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.Ixx.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Ixx.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.IxxIdx];
            end

            if self.paramDefs.Iyy.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.Iyy.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Iyy.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.IyyIdx];
            end

            if self.paramDefs.Izz.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.Izz.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Izz.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.IzzIdx];
            end

            if self.paramDefs.Ixy.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.Ixy.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Ixy.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.IxyIdx];
            end

            if self.paramDefs.Ixz.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.Ixz.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Ixz.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.IxzIdx];
            end

            if self.paramDefs.Iyz.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.Iyz.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Iyz.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.IyzIdx];
            end
            
            switch self.aeroModel
                case "constant"
                    % CX0
                    if self.paramDefs.CX0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CX0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CX0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CX0Idx];
                    end

                    % CX2
                    if self.paramDefs.CX2.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CX2.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CX2.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CX2Idx];
                    end

                    % CY0
                    if self.paramDefs.CY0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CY0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CY0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CY0Idx];
                    end

                    % CZ0
                    if self.paramDefs.CZ0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CZ0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CZ0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CZ0Idx];
                    end

                    % CNalpha0
                    if self.paramDefs.CNalpha0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CNalpha0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNalpha0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CNalpha0Idx];
                    end

                    % CNalpha2
                    if self.paramDefs.CNalpha2.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CNalpha2.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNalpha2.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CNalpha2Idx];
                    end

                    % CNpalpha0
                    if self.paramDefs.CNpalpha0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CNpalpha0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNpalpha0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CNpalpha0Idx];
                    end

                    % CNpalpha2
                    if self.paramDefs.CNpalpha2.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CNpalpha2.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNpalpha2.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CNpalpha2Idx];
                    end

                    % Cl0
                    if self.paramDefs.Cl0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.Cl0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cl0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.Cl0Idx];
                    end

                    % Clp
                    if self.paramDefs.Clp.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.Clp.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Clp.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.ClpIdx];
                    end

                    % Cldelta
                    if self.paramDefs.Cldelta.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.Cldelta.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cldelta.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CldeltaIdx];
                    end

                    % Cm0
                    if self.paramDefs.Cm0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.Cm0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cm0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.Cm0Idx];
                    end

                    % Cn0
                    if self.paramDefs.Cn0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.Cn0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cn0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.Cn0Idx];
                    end

                    % CMalpha0
                    if self.paramDefs.CMalpha0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CMalpha0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMalpha0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CMalpha0Idx];
                    end

                    % CMalpha2
                    if self.paramDefs.CMalpha2.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CMalpha2.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMalpha2.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CMalpha2Idx];
                    end

                    % CMpalpha0
                    if self.paramDefs.CMpalpha0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CMpalpha0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMpalpha0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CMpalpha0Idx];
                    end

                    % CMpalpha2
                    if self.paramDefs.CMpalpha2.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CMpalpha2.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMpalpha2.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CMpalpha2Idx];
                    end

                    % CMq
                    if self.paramDefs.CMq.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.CMq.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMq.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.CMqIdx];
                    end

                case "table"
                    % CX0
                    for i = 1:length(self.paramDefs.CX0.yValues)
                        if self.paramDefs.CX0.yIsEstimated(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CX0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CX0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CX0Table_CX0Idx + (i - 1)];
                        end
                    end

                    % CX2
                    for i = 1:length(self.paramDefs.CX2.yValues)
                        if self.paramDefs.CX2.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CX2.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CX2.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CX2Table_CX2Idx + (i - 1)];
                        end
                    end

                    % CY0
                    for i = 1:length(self.paramDefs.CY0.yValues)
                        if self.paramDefs.CY0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CY0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CY0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CY0Table_CY0Idx + (i - 1)];
                        end
                    end

                    % CZ0
                    for i = 1:length(self.paramDefs.CZ0.yValues)
                        if self.paramDefs.CZ0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CZ0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CZ0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CZ0Table_CZ0Idx + (i - 1)];
                        end
                    end

                    % CNalpha0
                    for i = 1:length(self.paramDefs.CNalpha0.yValues)
                        if self.paramDefs.CNalpha0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CNalpha0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNalpha0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CNalpha0Table_CNalpha0Idx + (i - 1)];
                        end
                    end

                    % CNalpha2
                    for i = 1:length(self.paramDefs.CNalpha2.yValues)
                        if self.paramDefs.CNalpha2.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CNalpha2.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNalpha2.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CNalpha2Table_CNalpha2Idx + (i - 1)];
                        end
                    end

                    % CNpalpha0
                    for i = 1:length(self.paramDefs.CNpalpha0.yValues)
                        if self.paramDefs.CNpalpha0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CNpalpha0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNpalpha0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CNpalpha0Table_CNpalpha0Idx + (i - 1)];
                        end
                    end

                    % CNpalpha2
                    for i = 1:length(self.paramDefs.CNpalpha2.yValues)
                        if self.paramDefs.CNpalpha2.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CNpalpha2.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CNpalpha2.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CNpalpha2Table_CNpalpha2Idx + (i - 1)];
                        end
                    end

                    % Cl0
                    for i = 1:length(self.paramDefs.Cl0.yValues)
                        if self.paramDefs.Cl0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.Cl0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cl0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.Cl0Table_Cl0Idx + (i - 1)];
                        end
                    end

                    % Clp
                    for i = 1:length(self.paramDefs.Clp.yValues)
                        if self.paramDefs.Clp.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.Clp.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Clp.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.ClpTable_ClpIdx + (i - 1)];
                        end
                    end

                    % Cldelta
                    for i = 1:length(self.paramDefs.Cldelta.yValues)
                        if self.paramDefs.Cldelta.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.Cldelta.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cldelta.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CldeltaTable_CldeltaIdx + (i - 1)];
                        end
                    end

                    % Cm0
                    for i = 1:length(self.paramDefs.Cm0.yValues)
                        if self.paramDefs.Cm0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.Cm0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cm0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.Cm0Table_Cm0Idx + (i - 1)];
                        end
                    end

                    % Cn0
                    for i = 1:length(self.paramDefs.Cn0.yValues)
                        if self.paramDefs.Cn0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.Cn0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.Cn0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.Cn0Table_Cn0Idx + (i - 1)];
                        end
                    end

                    % CMalpha0
                    for i = 1:length(self.paramDefs.CMalpha0.yValues)
                        if self.paramDefs.CMalpha0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CMalpha0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMalpha0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CMalpha0Table_CMalpha0Idx + (i - 1)];
                        end
                    end

                    % CMalpha2
                    for i = 1:length(self.paramDefs.CMalpha2.yValues)
                        if self.paramDefs.CMalpha2.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CMalpha2.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMalpha2.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CMalpha2Table_CMalpha2Idx + (i - 1)];
                        end
                    end

                    % CMpalpha0
                    for i = 1:length(self.paramDefs.CMpalpha0.yValues)
                        if self.paramDefs.CMpalpha0.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CMpalpha0.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMpalpha0.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CMpalpha0Table_CMpalpha0Idx + (i - 1)];
                        end
                    end

                    % CMpalpha2
                    for i = 1:length(self.paramDefs.CMpalpha2.yValues)
                        if self.paramDefs.CMpalpha2.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CMpalpha2.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMpalpha2.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CMpalpha2Table_CMpalpha2Idx + (i - 1)];
                        end
                    end

                    % CMq
                    for i = 1:length(self.paramDefs.CMq.yValues)
                        if self.paramDefs.CMq.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.CMq.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.CMq.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.CMqTable_CMqIdx + (i - 1)];
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
            self.paramDefs.Iyy       = ParamDef(props(6));
            self.paramDefs.Izz       = ParamDef(props(7));
            self.paramDefs.Ixy       = ParamDef(props(8));
            self.paramDefs.Ixz       = ParamDef(props(9));
            self.paramDefs.Iyz       = ParamDef(props(10));
        end


        function readAeroModelTablesFromFile(self, filePath)
            try
                aeroTable = readmatrix(filePath);
            catch
                error("Cannot read aero table CSV file. File either does not exist or is not formatted properly.")
            end

            machValues      = aeroTable(:, 1);
            CX0Values       = aeroTable(:, 2);
            CX2Values       = aeroTable(:, 3);
            CY0Values       = aeroTable(:, 4);
            CZ0Values       = aeroTable(:, 5);
            CNalpha0Values  = aeroTable(:, 6);
            CNalpha2Values  = aeroTable(:, 7);
            CNpalpha0Values = aeroTable(:, 8);
            CNpalpha2Values = aeroTable(:, 9);
            Cl0Values       = aeroTable(:, 10);
            ClpValues       = aeroTable(:, 11);
            CldeltaValues   = aeroTable(:, 12);
            Cm0Values       = aeroTable(:, 13);
            Cn0Values       = aeroTable(:, 14);
            CMalpha0Values  = aeroTable(:, 15);
            CMalpha2Values  = aeroTable(:, 16);
            CMpalpha0Values = aeroTable(:, 17);
            CMpalpha2Values = aeroTable(:, 18);
            CMqValues       = aeroTable(:, 19);

            self.paramDefs.CX0       = ParamTableDef(machValues, CX0Values);
            self.paramDefs.CX2       = ParamTableDef(machValues, CX2Values);
            self.paramDefs.CY0       = ParamTableDef(machValues, CY0Values);
            self.paramDefs.CZ0       = ParamTableDef(machValues, CZ0Values);
            self.paramDefs.CNalpha0  = ParamTableDef(machValues, CNalpha0Values);
            self.paramDefs.CNalpha2  = ParamTableDef(machValues, CNalpha2Values);
            self.paramDefs.CNpalpha0 = ParamTableDef(machValues, CNpalpha0Values);
            self.paramDefs.CNpalpha2 = ParamTableDef(machValues, CNpalpha2Values);
            self.paramDefs.Cl0       = ParamTableDef(machValues, Cl0Values);
            self.paramDefs.Clp       = ParamTableDef(machValues, ClpValues);
            self.paramDefs.Cldelta   = ParamTableDef(machValues, CldeltaValues);
            self.paramDefs.Cm0       = ParamTableDef(machValues, Cm0Values);
            self.paramDefs.Cn0       = ParamTableDef(machValues, Cn0Values);
            self.paramDefs.CMalpha0  = ParamTableDef(machValues, CMalpha0Values);
            self.paramDefs.CMalpha2  = ParamTableDef(machValues, CMalpha2Values);
            self.paramDefs.CMpalpha0 = ParamTableDef(machValues, CMpalpha0Values);
            self.paramDefs.CMpalpha2 = ParamTableDef(machValues, CMpalpha2Values);
            self.paramDefs.CMq       = ParamTableDef(machValues, CMqValues);
        end

        
        % Model methods ============================================================================

        function [CX0, CX2, CY0, CZ0, CNalpha0, CNalpha2, CNpalpha0, CNpalpha2, ...
                  Cl0, Clp, Cldelta, Cm0, Cn0, CMalpha0, CMalpha2, CMpalpha0, CMpalpha2, CMq] = constantAeroModel(self, ~)

            CX0 = self.params(self.CX0Idx);
            CX2 = self.params(self.CX2Idx);
            CY0 = self.params(self.CY0Idx);
            CZ0 = self.params(self.CZ0Idx);
            CNalpha0 = self.params(self.CNalpha0Idx);
            CNalpha2 = self.params(self.CNalpha2Idx);
            CNpalpha0 = self.params(self.CNpalpha0Idx);
            CNpalpha2 = self.params(self.CNpalpha2Idx);
            Cl0 = self.params(self.Cl0Idx);
            Clp = self.params(self.ClpIdx);
            Cldelta = self.params(self.CldeltaIdx);
            Cm0 = self.params(self.Cm0Idx);
            Cn0 = self.params(self.Cn0Idx);
            CMalpha0 = self.params(self.CMalpha0Idx);
            CMalpha2 = self.params(self.CMalpha2Idx);
            CMpalpha0 = self.params(self.CMpalpha0Idx);
            CMpalpha2 = self.params(self.CMpalpha2Idx);
            CMq = self.params(self.CMqIdx);
        end


        function [CX0, CX2, CY0, CZ0, CNalpha0, CNalpha2, CNpalpha0, CNpalpha2, ...
                  Cl0, Clp, Cldelta, Cm0, Cn0, CMalpha0, CMalpha2, CMpalpha0, CMpalpha2, CMq] = tableAeroModel(self, mach)
            
            minMach = self.params(self.CX0Table_Mach0Idx);

            if mach < minMach
                CX0       = self.params(self.CX0Table_CX0Idx);
                CX2       = self.params(self.CX2Table_CX2Idx);
                CY0       = self.params(self.CY0Table_CY0Idx);
                CZ0       = self.params(self.CZ0Table_CZ0Idx);
                CNalpha0  = self.params(self.CNalpha0Table_CNalpha0Idx);
                CNalpha2  = self.params(self.CNalpha2Table_CNalpha2Idx);
                CNpalpha0 = self.params(self.CNpalpha0Table_CNpalpha0Idx);
                CNpalpha2 = self.params(self.CNpalpha2Table_CNpalpha2Idx);
                Cl0       = self.params(self.Cl0Table_Cl0Idx);
                Clp       = self.params(self.ClpTable_ClpIdx);
                Cldelta   = self.params(self.CldeltaTable_CldeltaIdx);
                Cm0       = self.params(self.Cm0Table_Cm0Idx);
                Cn0       = self.params(self.Cn0Table_Cn0Idx);
                CMalpha0  = self.params(self.CMalpha0Table_CMalpha0Idx);
                CMalpha2  = self.params(self.CMalpha2Table_CMalpha2Idx);
                CMpalpha0 = self.params(self.CMpalpha0Table_CMpalpha0Idx);
                CMpalpha2 = self.params(self.CMpalpha2Table_CMpalpha2Idx);
                CMq       = self.params(self.CMqTable_CMqIdx);

            else
                dmach = self.params(self.CX0Table_Mach0Idx + 1) - self.params(self.CX0Table_Mach0Idx);  % TODO: Non-uniform table
    
                CX0 = 0;
                for i = 0:(self.CX0Table_Len - 1)
                    mach_i = self.params(self.CX0Table_Mach0Idx + i);
                    CX0_i = self.params(self.CX0Table_CX0Idx + i);
    
                    CX0 = CX0 + CX0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CX2 = 0;
                for i = 0:(self.CX2Table_Len - 1)
                    mach_i = self.params(self.CX2Table_Mach0Idx + i);
                    CX2_i = self.params(self.CX2Table_CX2Idx + i);
    
                    CX2 = CX2 + CX2_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CY0 = 0;
                for i = 0:(self.CY0Table_Len - 1)
                    mach_i = self.params(self.CY0Table_Mach0Idx + i);
                    CY0_i = self.params(self.CY0Table_CY0Idx + i);
    
                    CY0 = CY0 + CY0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CZ0 = 0;
                for i = 0:(self.CZ0Table_Len - 1)
                    mach_i = self.params(self.CZ0Table_Mach0Idx + i);
                    CZ0_i = self.params(self.CZ0Table_CZ0Idx + i);
    
                    CZ0 = CZ0 + CZ0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CNalpha0 = 0;
                for i = 0:(self.CNalpha0Table_Len - 1)
                    mach_i = self.params(self.CNalpha0Table_Mach0Idx + i);
                    CNalpha0_i = self.params(self.CNalpha0Table_CNalpha0Idx + i);
    
                    CNalpha0 = CNalpha0 + CNalpha0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CNalpha2 = 0;
                for i = 0:(self.CNalpha2Table_Len - 1)
                    mach_i = self.params(self.CNalpha2Table_Mach0Idx + i);
                    CNalpha2_i = self.params(self.CNalpha2Table_CNalpha2Idx + i);
    
                    CNalpha2 = CNalpha2 + CNalpha2_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CNpalpha0 = 0;
                for i = 0:(self.CNpalpha0Table_Len - 1)
                    mach_i = self.params(self.CNpalpha0Table_Mach0Idx + i);
                    CNpalpha0_i = self.params(self.CNpalpha0Table_CNpalpha0Idx + i);
    
                    CNpalpha0 = CNpalpha0 + CNpalpha0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CNpalpha2 = 0;
                for i = 0:(self.CNpalpha2Table_Len - 1)
                    mach_i = self.params(self.CNpalpha2Table_Mach0Idx + i);
                    CNpalpha2_i = self.params(self.CNpalpha2Table_CNpalpha2Idx + i);
    
                    CNpalpha2 = CNpalpha2 + CNpalpha2_i * self.linearKernel((mach - mach_i) / dmach);
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
    
                Cm0 = 0;
                for i = 0:(self.Cm0Table_Len - 1)
                    mach_i = self.params(self.Cm0Table_Mach0Idx + i);
                    Cm0_i = self.params(self.Cm0Table_Cm0Idx + i);
    
                    Cm0 = Cm0 + Cm0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                Cn0 = 0;
                for i = 0:(self.Cn0Table_Len - 1)
                    mach_i = self.params(self.Cn0Table_Mach0Idx + i);
                    Cn0_i = self.params(self.Cn0Table_Cn0Idx + i);
    
                    Cn0 = Cn0 + Cn0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CMalpha0 = 0;
                for i = 0:(self.CMalpha0Table_Len - 1)
                    mach_i = self.params(self.CMalpha0Table_Mach0Idx + i);
                    CMalpha0_i = self.params(self.CMalpha0Table_CMalpha0Idx + i);
    
                    CMalpha0 = CMalpha0 + CMalpha0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CMalpha2 = 0;
                for i = 0:(self.CMalpha2Table_Len - 1)
                    mach_i = self.params(self.CMalpha2Table_Mach0Idx + i);
                    CMalpha2_i = self.params(self.CMalpha2Table_CMalpha2Idx + i);
    
                    CMalpha2 = CMalpha2 + CMalpha2_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CMpalpha0 = 0;
                for i = 0:(self.CMpalpha0Table_Len - 1)
                    mach_i = self.params(self.CMpalpha0Table_Mach0Idx + i);
                    CMpalpha0_i = self.params(self.CMpalpha0Table_CMpalpha0Idx + i);
    
                    CMpalpha0 = CMpalpha0 + CMpalpha0_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CMpalpha2 = 0;
                for i = 0:(self.CMpalpha2Table_Len - 1)
                    mach_i = self.params(self.CMpalpha2Table_Mach0Idx + i);
                    CMpalpha2_i = self.params(self.CMpalpha2Table_CMpalpha2Idx + i);
    
                    CMpalpha2 = CMpalpha2 + CMpalpha2_i * self.linearKernel((mach - mach_i) / dmach);
                end
    
                CMq = 0;
                for i = 0:(self.CMqTable_Len - 1)
                    mach_i = self.params(self.CMqTable_Mach0Idx + i);
                    CMq_i = self.params(self.CMqTable_CMqIdx + i);
    
                    CMq = CMq + CMq_i * self.linearKernel((mach - mach_i) / dmach);
                end
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

        function set.consideredParams(self, consideredParams)
            if Settings.VALIDATE_FLAG
                self.consideredParams = Validator.validateType(consideredParams, "double");
            else
                self.consideredParams = consideredParams;
            end

            self.nConsideredParams = length(consideredParams);
        end

        function set.consideredParamCovar(self, consideredParamCovar)
            if Settings.VALIDATE_FLAG
                consideredParamCovar = Validator.validateType(consideredParamCovar, "double");
                self.consideredParamCovar = Validator.validateSize(consideredParamCovar, [self.nConsideredParams, self.nConsideredParams]);
            else
                self.consideredParamCovar = consideredParamCovar;
            end
        end

        function set.consideredParamIdxs(self, consideredParamIdxs)
            if Settings.VALIDATE_FLAG
                self.consideredParamIdxs = Validator.validateType(consideredParamIdxs, "double");
            else
                self.consideredParamIdxs = consideredParamIdxs;
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
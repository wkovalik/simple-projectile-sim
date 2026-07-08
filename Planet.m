classdef Planet < handle
    % TODO: Make a static Kernel class. Then set kernel property to choose which kernel to use
    
    properties
        paramDefs

        params = [];

        estimatedParams = [];
        estimatedParamCovar = [];

        consideredParams = [];
        consideredParamCovar = [];

        gravityModel
        atmosphereModel
        windModel
        windModelKernel
    end

    properties (SetAccess = protected)
        nParams = 0;

        nEstimatedParams = 0;
        estimatedParamIdxs = [];

        nConsideredParams = 0;
        consideredParamIdxs = [];

        gIdx = 0;
        
        rhoIdx = 0;
        rho0Idx = 0;
        rhoTable_h0Idx = 0;
        rhoTable_rho0Idx = 0;
        rhoTable_Len = 0;

        aIdx = 0;
        
        HIdx = 0;
        
        vWindxIdx = 0;
        vWindxTable_h0Idx = 0;
        vWindxTable_vWindx0Idx = 0;
        vWindxTable_Len = 0;

        vWindyIdx = 0;
        vWindyTable_h0Idx = 0;
        vWindyTable_vWindy0Idx = 0;
        vWindyTable_Len = 0;

        isGravityModelInit = false;
        isAtmosphereModelInit = false;
        isWindModelInit = false;

        computeGravity
        computeAtmosphere
        computeWind

        windModelKernelID
    end

    properties (Constant)
        VALID_GRAVITY_MODELS = ["constant"];
        VALID_ATMOSPHERE_MODELS = ["constant", "exponential", "table"];
        VALID_WIND_MODELS = ["constant", "table"];
        VALID_WIND_MODEL_KERNELS = ["constant", "linear"];

        DEFAULT_WIND_MODEL_KERNEL = "linear";
    end

    properties (Constant, Abstract)
        DEFAULT_G
        
        DEFAULT_RHO
        DEFAULT_RHO0
        DEFAULT_RHO_TABLE_X
        DEFAULT_RHO_TABLE_Y

        % DEFAULT_R
        % DEFAULT_T
        DEFAULT_A  % TODO: Compute using R and T as params instead
        
        DEFAULT_H

        DEFAULT_VWINDX
        DEFAULT_VWINDX_TABLE_X
        DEFAULT_VWINDX_TABLE_Y
        DEFAULT_VWINDY
        DEFAULT_VWINDY_TABLE_X
        DEFAULT_VWINDY_TABLE_Y
    end


    methods
        % Constructor ==============================================================================

        function self = Planet(gravityModel, atmosphereModel, windModel)
            if nargin == 0
                self.gravityModel = "constant";
                self.atmosphereModel = "constant";
                self.windModel = "constant";

            elseif nargin == 1
                self.gravityModel = gravityModel;
                self.atmosphereModel = "constant";
                self.windModel = "constant";
            
            elseif nargin == 2
                self.gravityModel = gravityModel;
                self.atmosphereModel = atmosphereModel;
                self.windModel = "constant";

            elseif nargin == 3
                self.gravityModel = gravityModel;
                self.atmosphereModel = atmosphereModel;
                self.windModel = windModel;

            else
                error("Too many input arguments.")
            end
            
            self.update();
        end


        % Update methods ===========================================================================

        function update(self)
            self.updateModels();
            self.updateParams();
            self.updateEstimatedParams();
            self.updateConsideredParams();
        end


        function updateModels(self)
            self.updateGravityModel();
            self.updateAtmosphereModel();
            self.updateWindModel();
        end


        function updateGravityModel(self)
            switch self.gravityModel
                case "constant"
                    self.computeGravity = @self.constantGravityModel;
                    
                    if ~self.isGravityModelInit
                        self.paramDefs.g = ParamDef(self.DEFAULT_G);
                    end

                otherwise
                    error("Invalid gravity model: %s.", self.gravityModel)
            end

            self.isGravityModelInit = true;
        end


        function updateAtmosphereModel(self)
            switch self.atmosphereModel
                case "constant"
                    self.computeAtmosphere = @self.constantAtmosphereModel;
                    
                    if ~self.isAtmosphereModelInit
                        self.paramDefs.rho = ParamDef(self.DEFAULT_RHO);
                        self.paramDefs.a = ParamDef(self.DEFAULT_A);
                    end

                case "exponential"
                    self.computeAtmosphere = @self.exponentialAtmosphereModel;
                    
                    if ~self.isAtmosphereModelInit
                        self.paramDefs.rho0 = ParamDef(self.DEFAULT_RHO0);
                        self.paramDefs.H = ParamDef(self.DEFAULT_H);
                        self.paramDefs.a = ParamDef(self.DEFAULT_A);
                    end

                case "table"
                    self.computeAtmosphere = @self.tableAtmosphereModel;
                    
                    if ~self.isAtmosphereModelInit
                        self.paramDefs.rho = ParamTableDef(self.DEFAULT_RHO_TABLE_X, self.DEFAULT_RHO_TABLE_Y);
                        self.paramDefs.a = ParamDef(self.DEFAULT_A);
                    end

                otherwise
                    error("Invalid atmosphere model: %s.", self.atmosphereModel)
            end

            self.isAtmosphereModelInit = true;
        end


        function updateWindModel(self)
            switch self.windModel
                case "constant"
                    self.computeWind = @self.constantWindModel;
                    
                    if ~self.isWindModelInit
                        self.paramDefs.vWindx = ParamDef(self.DEFAULT_VWINDX);
                        self.paramDefs.vWindy = ParamDef(self.DEFAULT_VWINDY);
                    end

                case "table"
                    self.computeWind = @self.tableWindModel;
                    
                    if ~self.isWindModelInit
                        self.paramDefs.vWindx = ParamTableDef(self.DEFAULT_VWINDX_TABLE_X, self.DEFAULT_VWINDX_TABLE_Y);
                        self.paramDefs.vWindy = ParamTableDef(self.DEFAULT_VWINDY_TABLE_X, self.DEFAULT_VWINDY_TABLE_Y);

                        self.windModelKernel = self.DEFAULT_WIND_MODEL_KERNEL;
                        self.updateWindModelKernelID();
                    end

                otherwise
                    error("Invalid wind model: %s.", self.windModel)
            end
            
            self.isWindModelInit = true;
        end


        function updateWindModelKernelID(self)
            switch self.windModelKernel
                case "constant"
                    self.windModelKernelID = 1;

                case "linear"
                    self.windModelKernelID = 2;
                    
                otherwise
                    error("Invalid wind model kernel: %s.", self.windModelKernel)
            end
        end


        function updateParams(self)
            % See Note 1 regarding parameter indexing

            self.params = [];
            
            self.gIdx = 0;
        
            self.rhoIdx = 0;
            self.rho0Idx = 0;
            self.rhoTable_h0Idx = 0;
            self.rhoTable_rho0Idx = 0;
            self.rhoTable_Len = 0;
            
            self.HIdx = 0;

            self.aIdx = 0;
            
            self.vWindxIdx = 0;
            self.vWindxTable_h0Idx = 0;
            self.vWindxTable_vWindx0Idx = 0;
            self.vWindxTable_Len = 0;
    
            self.vWindyIdx = 0;
            self.vWindyTable_h0Idx = 0;
            self.vWindyTable_vWindy0Idx = 0;
            self.vWindyTable_Len = 0;

            switch self.gravityModel
                case "constant"
                    % g
                    self.gIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.g.value];
                
                otherwise
                    error("Invalid gravity model: %s.", self.gravityModel)
            end

            switch self.atmosphereModel
                case "constant"
                    % rho
                    self.rhoIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.rho.value];
                    
                    % a
                    self.aIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.a.value];

                case "exponential"
                    % rho0
                    self.rho0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.rho0.value];
                    
                    % H
                    self.HIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.H.value];
                    
                    % a
                    self.aIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.a.value];

                case "table"
                    % rho
                    self.rhoTable_h0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.rho.xValues];
                    
                    self.rhoTable_rho0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.rho.yValues];

                    self.rhoTable_Len = self.paramDefs.rho.nValues;

                    % a
                    self.aIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.a.value];

                otherwise
                    error("Invalid atmosphere model: %s.", self.atmosphereModel)
            end

            switch self.windModel
                case "constant"
                    % vWindx
                    self.vWindxIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.vWindx.value];
                    
                    % vWindy
                    self.vWindyIdx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.vWindy.value];

                case "table"
                    % vWindx
                    self.vWindxTable_h0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.vWindx.xValues];
                    
                    self.vWindxTable_vWindx0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.vWindx.yValues];

                    self.vWindxTable_Len = self.paramDefs.vWindy.nValues;
                    
                    % vWindy
                    self.vWindyTable_h0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.vWindy.xValues];
                    
                    self.vWindyTable_vWindy0Idx = self.nParams + 1;
                    self.params = [self.params; self.paramDefs.vWindy.yValues];

                    self.vWindyTable_Len = self.paramDefs.vWindy.nValues;
                
                otherwise
                    error("Invalid wind model: %s.", self.windModel)
            end
        end


        function updateEstimatedParams(self)
            self.estimatedParams = [];
            self.estimatedParamCovar = [];
            self.estimatedParamIdxs = [];

            switch self.gravityModel
                % g
                case "constant"
                    if self.paramDefs.g.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.g.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.g.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.gIdx];
                    end
                
                otherwise
                    error("Invalid gravity model: %s.", self.gravityModel)
            end

            switch self.atmosphereModel
                case "constant"
                    % rho
                    if self.paramDefs.rho.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.rho.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.rho.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.rhoIdx];
                    end
                    
                    % a
                    if self.paramDefs.a.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.a.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.a.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.aIdx];
                    end
                
                case "exponential"
                    % rho0
                    if self.paramDefs.rho0.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.rho0.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.rho0.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.rho0Idx];
                    end
                    
                    % H
                    if self.paramDefs.H.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.H.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.H.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.HIdx];
                    end
                    
                    % a
                    if self.paramDefs.a.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.a.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.a.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.aIdx];
                    end

                case "table"
                    % rho
                    for i = 1:length(self.paramDefs.rho.yValues)
                        if self.paramDefs.rho.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.rho.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.rho.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.rhoTable_rho0Idx + (i - 1)];
                        end
                    end

                    % a
                    if self.paramDefs.a.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.a.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.a.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.aIdx];
                    end
                
                otherwise
                    error("Invalid atmosphere model: %s.", self.atmosphereModel)
            end

            switch self.windModel
                case "constant"
                    % vWindx
                    if self.paramDefs.vWindx.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.vWindx.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.vWindx.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.vWindxIdx];
                    end
                    
                    % vWindy
                    if self.paramDefs.vWindy.isEstimated
                        self.estimatedParams = [self.estimatedParams; self.paramDefs.vWindy.value];
                        self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.vWindy.covar);
                        self.estimatedParamIdxs = [self.estimatedParamIdxs; self.vWindyIdx];
                    end

                case "table"
                    % vWindx
                    for i = 1:length(self.paramDefs.vWindx.yValues)
                        if self.paramDefs.vWindx.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.vWindx.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.vWindx.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.vWindxTable_vWindx0Idx + (i - 1)];
                        end
                    end
                    
                    % vWindy
                    for i = 1:length(self.paramDefs.vWindy.yValues)
                        if self.paramDefs.vWindy.yIsEstimated(i)
                            self.estimatedParams = [self.estimatedParams; self.paramDefs.vWindy.yValues(i)];
                            self.estimatedParamCovar = blkdiag(self.estimatedParamCovar, self.paramDefs.vWindy.yCovars(i));
                            self.estimatedParamIdxs = [self.estimatedParamIdxs; self.vWindyTable_vWindy0Idx + (i - 1)];
                        end
                    end

                otherwise
                    error("Invalid wind model: %s.", self.windModel)
            end
        end


        function updateConsideredParams(self)
            self.consideredParams = [];
            self.consideredParamCovar = [];
            self.consideredParamIdxs = [];

            switch self.gravityModel
                % g
                case "constant"
                    if self.paramDefs.g.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.g.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.g.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.gIdx];
                    end
                
                otherwise
                    error("Invalid gravity model: %s.", self.gravityModel)
            end

            switch self.atmosphereModel
                case "constant"
                    % rho
                    if self.paramDefs.rho.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.rho.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.rho.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.rhoIdx];
                    end
                    
                    % a
                    if self.paramDefs.a.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.a.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.a.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.aIdx];
                    end
                
                case "exponential"
                    % rho0
                    if self.paramDefs.rho0.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.rho0.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.rho0.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.rho0Idx];
                    end
                    
                    % H
                    if self.paramDefs.H.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.H.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.H.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.HIdx];
                    end
                    
                    % a
                    if self.paramDefs.a.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.a.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.a.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.aIdx];
                    end

                case "table"
                    % rho
                    for i = 1:length(self.paramDefs.rho.yValues)
                        if self.paramDefs.rho.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.rho.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.rho.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.rhoTable_rho0Idx + (i - 1)];
                        end
                    end

                    % a
                    if self.paramDefs.a.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.a.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.a.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.aIdx];
                    end
                
                otherwise
                    error("Invalid atmosphere model: %s.", self.atmosphereModel)
            end

            switch self.windModel
                case "constant"
                    % vWindx
                    if self.paramDefs.vWindx.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.vWindx.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.vWindx.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.vWindxIdx];
                    end
                    
                    % vWindy
                    if self.paramDefs.vWindy.isConsidered
                        self.consideredParams = [self.consideredParams; self.paramDefs.vWindy.value];
                        self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.vWindy.covar);
                        self.consideredParamIdxs = [self.consideredParamIdxs; self.vWindyIdx];
                    end

                case "table"
                    % vWindx
                    for i = 1:length(self.paramDefs.vWindx.yValues)
                        if self.paramDefs.vWindx.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.vWindx.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.vWindx.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.vWindxTable_vWindx0Idx + (i - 1)];
                        end
                    end
                    
                    % vWindy
                    for i = 1:length(self.paramDefs.vWindy.yValues)
                        if self.paramDefs.vWindy.yIsConsidered(i)
                            self.consideredParams = [self.consideredParams; self.paramDefs.vWindy.yValues(i)];
                            self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.vWindy.yCovars(i));
                            self.consideredParamIdxs = [self.consideredParamIdxs; self.vWindyTable_vWindy0Idx + (i - 1)];
                        end
                    end

                otherwise
                    error("Invalid wind model: %s.", self.windModel)
            end
        end


        function readDensityTableFromFile(self, filePath, usePerturbedValues)
            if nargin < 3
                usePerturbedValues = false;
            end

            try
                atmTable = readmatrix(filePath);
            catch
                error("Cannot read atmosphere table CSV file. File either does not exist or is not formatted properly.")
            end

            hValues = atmTable(:, 2);
            hValues = hValues * 1E+03;  % (km) to (m)

            if ~usePerturbedValues
                rhoValues    =  atmTable(:, 10);
            else
                rhoValues    =  atmTable(:, 23);
            end

            self.paramDefs.rho    = ParamTableDef(hValues, rhoValues);
        end


        function readWindTablesFromFile(self, filePath, usePerturbedValues)
            if nargin < 3
                usePerturbedValues = false;
            end

            try
                atmTable = readmatrix(filePath);
            catch
                error("Cannot read atmosphere table CSV file. File either does not exist or is not formatted properly.")
            end

            hValues = atmTable(:, 2);
            hValues = hValues * 1E+03;  % (km) to (m)

            if ~usePerturbedValues
                vWindxValues =  atmTable(:, 32);   % Assume x-axis aligned due East (so East wind is positive)
                vWindyValues = -atmTable(:, 33);   % Assume y-axis aligned due South (so North wind is negative)
            else
                vWindxValues =  atmTable(:, 38);
                vWindyValues = -atmTable(:, 39);
            end

            self.paramDefs.vWindx = ParamTableDef(hValues, vWindxValues);
            self.paramDefs.vWindy = ParamTableDef(hValues, vWindyValues);
        end

        
        % Model methods ============================================================================

        function g = constantGravityModel(self)
            g = self.params(self.gIdx);
        end

        
        function [rho, a] = constantAtmosphereModel(self, ~)
            rho = self.params(self.rhoIdx);
            a = self.params(self.aIdx);
        end


        function [rho, a] = exponentialAtmosphereModel(self, h)
            rho0 = self.params(self.rho0Idx);
            H = self.params(self.HIdx);
            a = self.params(self.aIdx);
            
            rho = rho0 * exp(-h / H);
        end

        function [rho, a] = tableAtmosphereModel(self, h)
            a = self.params(self.aIdx);

            % TODO: Non-uniform table
            dh = self.params(self.rhoTable_h0Idx + 1) - self.params(self.rhoTable_h0Idx);

            rho = 0;
            for i = 0:(self.rhoTable_Len - 1)
                h_i = self.params(self.rhoTable_h0Idx + i);
                normalized_dh_i = (h - h_i) / dh;

                rho_i = self.params(self.rhoTable_rho0Idx + i);

                rho = rho + rho_i * self.linearKernel(normalized_dh_i);
            end
        end


        function [vWindx, vWindy] = constantWindModel(self, ~)
            vWindx = self.params(self.vWindxIdx);
            vWindy = self.params(self.vWindyIdx);
        end


        function [vWindx, vWindy] = tableWindModel(self, h)
            % TODO: Non-uniform table
            dh = self.params(self.vWindxTable_h0Idx + 1) - self.params(self.vWindxTable_h0Idx);  % Using vWindx table since entire table uses same height points

            vWindx = 0;
            vWindy = 0;

            for i = 0:(self.vWindxTable_Len - 1)  % Using vWindx table since entire table uses same height points
                h_i = self.params(self.vWindxTable_h0Idx + i);
                normalized_dh_i = (h - h_i) / dh;

                vWindx_i = self.params(self.vWindxTable_vWindx0Idx + i);
                vWindy_i = self.params(self.vWindyTable_vWindy0Idx + i);

                if self.windModelKernelID == 1
                    vWindx = vWindx + vWindx_i * self.constantKernel(normalized_dh_i);
                    vWindy = vWindy + vWindy_i * self.constantKernel(normalized_dh_i);
                
                elseif self.windModelKernelID == 2
                    vWindx = vWindx + vWindx_i * self.linearKernel(normalized_dh_i);
                    vWindy = vWindy + vWindy_i * self.linearKernel(normalized_dh_i);

                end
            end
        end


        function k = constantKernel(~, x)
            if (-0.5 <= x) && (x < 0.5)
                k = 1;
            else
                k = 0;
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

        function set.paramDefs(self, paramDefs)
            if Settings.VALIDATE_FLAG
                self.paramDefs = Validator.validateFieldTypes(paramDefs, ["ParamDef", "ParamTableDef"]);
            else
                self.paramDefs = paramDefs;
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

        function set.gravityModel(self, gravityModel)
            if Settings.VALIDATE_FLAG
                self.gravityModel = Validator.validateString(gravityModel, self.VALID_GRAVITY_MODELS);
            else
                self.gravityModel = gravityModel;
            end

            self.isGravityModelInit = false;
            self.updateGravityModel();
        end

        function set.atmosphereModel(self, atmosphereModel)
            if Settings.VALIDATE_FLAG
                self.atmosphereModel = Validator.validateString(atmosphereModel, self.VALID_ATMOSPHERE_MODELS);
            else
                self.atmosphereModel = atmosphereModel;
            end

            self.isAtmosphereModelInit = false;
            self.updateAtmosphereModel();
        end

        function set.windModel(self, windModel)
            if Settings.VALIDATE_FLAG
                self.windModel = Validator.validateString(windModel, self.VALID_WIND_MODELS);
            else
                self.windModel = windModel;
            end

            self.isWindModelInit = false;
            self.updateWindModel();
        end

        function set.windModelKernel(self, windModelKernel)
            if Settings.VALIDATE_FLAG
                self.windModelKernel = Validator.validateString(windModelKernel, self.VALID_WIND_MODEL_KERNELS);
            else
                self.windModelKernel = windModelKernel;
            end

            self.updateWindModelKernelID();
        end

        function set.computeGravity(self, gravityModelFn)
            if Settings.VALIDATE_FLAG
                self.computeGravity = Validator.validateType(gravityModelFn, "function_handle");
            else
                self.computeGravity = gravityModelFn;
            end
        end

        function set.computeAtmosphere(self, atmosphereModelFn)
            if Settings.VALIDATE_FLAG
                self.computeAtmosphere = Validator.validateType(atmosphereModelFn, "function_handle");
            else
                self.computeAtmosphere = atmosphereModelFn;
            end
        end

        function set.computeWind(self, windModelFn)
            if Settings.VALIDATE_FLAG
                self.computeWind = Validator.validateType(windModelFn, "function_handle");
            else
                self.computeWind = windModelFn;
            end
        end
    end
end


% Note 1
%
% Parameter values are obtained by directly getting the value from the self.params vector using the
% respective parameter index. For example:
% 
% | g = self.params(self.gIdx);
% 
% This is extremely fast. A possible alternative is to abstract this direct indexing call behind a
% dependent "parameter" with a getter. For example, define:
%
% | properties (Dependent)
% |     g
% | end
% 
% | methods
% |     function g = get.g(self)
% |         g = self.params(self.gIdx);
% |     end
% | end
%
% Then, to get the parameter value:
%
% | g = self.g
%
% which calls the getter. This results is arguably more readable code, at the cost of a significant
% performance hit (~15% slower).
%
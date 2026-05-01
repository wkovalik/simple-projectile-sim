classdef Planet < handle
    % TODO: Make a static Kernel class. Then set kernel property to choose which kernel to use
    
    properties
        paramDefs

        params = [];

        estimatedParams = [];
        estimatedParamCovar = [];

        gravityModel
        atmosphereModel
        windModel
    end

    properties (SetAccess = private)
        nParams = 0;

        nEstimatedParams = 0;
        estimatedParamIdxs = [];

        gIdx = 0;
        
        rhoIdx = 0;
        aIdx = 0;
        rho0Idx = 0;
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
    end

    properties (Constant)
        VALID_GRAVITY_MODELS = ["constant"];
        VALID_ATMOSPHERE_MODELS = ["constant", "exponential"];  % TODO: "table"
        VALID_WIND_MODELS = ["constant", "table"];
    end

    properties (Constant, Abstract)
        DEFAULT_G
        
        DEFAULT_RHO
        DEFAULT_RHO0
        DEFAULT_H
        % DEFAULT_R
        % DEFAULT_T
        DEFAULT_A  % TODO: Compute using R and T as params instead

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
                    end

                otherwise
                    error("Invalid wind model: %s.", self.windModel)
            end
            
            self.isWindModelInit = true;
        end


        function updateParams(self)
            % See Note 1 regarding parameter indexing

            self.params = [];
            
            self.gIdx = 0;
        
            self.rhoIdx = 0;
            self.aIdx = 0;
            self.rho0Idx = 0;
            self.HIdx = 0;
            
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


        function [vWindx, vWindy] = constantWindModel(self, ~)
            vWindx = self.params(self.vWindxIdx);
            vWindy = self.params(self.vWindyIdx);
        end


        function [vWindx, vWindy] = tableWindModel(self, h)
            dh = self.params(self.vWindxTable_h0Idx + 1) - self.params(self.vWindxTable_h0Idx);

            vWindx = 0;
            for i = 0:(self.vWindxTable_Len - 1)
                h_i = self.params(self.vWindxTable_h0Idx + i);
                vWindx_i = self.params(self.vWindxTable_vWindx0Idx + i);

                vWindx = vWindx + vWindx_i * self.linearKernel(h - h_i, dh);
            end

            dh = self.params(self.vWindyTable_h0Idx + 1) - self.params(self.vWindyTable_h0Idx);

            vWindy = 0;
            for i = 0:(self.vWindxTable_Len - 1)
                h_i = self.params(self.vWindyTable_h0Idx + i);
                vWindy_i = self.params(self.vWindyTable_vWindy0Idx + i);

                vWindy = vWindy + vWindy_i * self.linearKernel(h - h_i, dh);
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
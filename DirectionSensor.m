classdef DirectionSensor < Sensor
    properties (SetAccess = private)
        xIdx = 0;
        yIdx = 0;
        zIdx = 0;
        azBiasIdx = 0;
        elBiasIdx = 0;
    end

    properties (Constant)
        nMeas = 2;

        DEFAULT_X = 0;
        DEFAULT_Y = 0;
        DEFAULT_Z = 0;
        DEFAULT_AZBIAS = 0;
        DEFAULT_ELBIAS = 0;
    end

    
    methods
        % Constructor ==============================================================================

        function self = DirectionSensor()
            self = self@Sensor();

            self.paramDefs.x = ParamDef(self.DEFAULT_X);
            self.paramDefs.y = ParamDef(self.DEFAULT_Y);
            self.paramDefs.z = ParamDef(self.DEFAULT_Z);
            self.paramDefs.azBias = ParamDef(self.DEFAULT_AZBIAS);
            self.paramDefs.elBias = ParamDef(self.DEFAULT_ELBIAS);

            self.updateParams();
        end

        
        % Update methods ===========================================================================

        function updateParams(self)
            self.params = [];

            self.xIdx = 0;
            self.yIdx = 0;
            self.zIdx = 0;
            self.azBiasIdx = 0;
            self.elBiasIdx = 0;

            self.xIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.x.value];

            self.yIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.y.value];

            self.zIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.z.value];

            self.azBiasIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.azBias.value];

            self.elBiasIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.elBias.value];
        end


        function updateConsideredParams(self)
            self.consideredParams = [];
            self.consideredParamCovar = [];
            self.consideredParamIdxs = [];

            % x
            if self.paramDefs.x.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.x.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.x.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.xIdx];
            end

            % y
            if self.paramDefs.y.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.y.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.y.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.yIdx];
            end

            % z
            if self.paramDefs.z.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.z.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.z.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.zIdx];
            end

            % azBias
            if self.paramDefs.azBias.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.azBias.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.azBias.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.azBiasIdx];
            end

            % elBias
            if self.paramDefs.elBias.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.elBias.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.elBias.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.elBiasIdx];
            end
        end
        

        % Measurement methods ======================================================================

        function dir = takeMeasurement(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params(self.xIdx);
            ySensor = self.params(self.yIdx);
            zSensor = self.params(self.zIdx);
            azBias = self.params(self.azBiasIdx);
            elBias = self.params(self.elBiasIdx);
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;

            % Compute true relative downrange
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            % Compute true azimuth and elevation (with bias)
            az = atan2(dySensor, dxSensor) + azBias;
            el = -atan2(dzSensor, drSensor) + elBias;
            
            dir = [az; el];
            
            % Compute noisy azimuth and elevation measurements
            epsDir = self.measNoiseCovarSqrt * randn(self.nMeas, 1);

            dir = dir + epsDir;
        end


        function dir = computeMeasurement(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params(self.xIdx);
            ySensor = self.params(self.yIdx);
            zSensor = self.params(self.zIdx);
            azBias = self.params(self.azBiasIdx);
            elBias = self.params(self.elBiasIdx);
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;

            % Compute relative downrange
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            % Compute azimuth and elevation (with bias)
            az = atan2(dySensor, dxSensor) + azBias;
            el = -atan2(dzSensor, drSensor) + elBias;
            
            dir = [az; el];
        end
    end
end
classdef RangeSensor < Sensor
    properties (SetAccess = private)
        xIdx = 0;
        yIdx = 0;
        zIdx = 0;
        biasIdx = 0;
    end

    properties (Constant)
        nMeas = 1;

        DEFAULT_X = 0;
        DEFAULT_Y = 0;
        DEFAULT_Z = 0;
        DEFAULT_BIAS = 0;
    end


    methods
        % Constructor ==============================================================================

        function self = RangeSensor()
            self = self@Sensor();

            self.paramDefs.x = ParamDef(self.DEFAULT_X);
            self.paramDefs.y = ParamDef(self.DEFAULT_Y);
            self.paramDefs.z = ParamDef(self.DEFAULT_Z);
            self.paramDefs.bias = ParamDef(self.DEFAULT_BIAS);

            self.updateParams();
        end

        
        % Update methods ===========================================================================

        function updateParams(self)
            self.params = [];

            self.xIdx = 0;
            self.yIdx = 0;
            self.zIdx = 0;
            self.biasIdx = 0;

            self.xIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.x.value];

            self.yIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.y.value];

            self.zIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.z.value];

            self.biasIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.bias.value];
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

            % bias
            if self.paramDefs.bias.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.bias.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.bias.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.biasIdx];
            end
        end
        

        % Measurement methods ======================================================================
        
        function R = takeMeasurement(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params(self.xIdx);
            ySensor = self.params(self.yIdx);
            zSensor = self.params(self.zIdx);
            bias = self.params(self.biasIdx);
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;
            
            % Compute true range (with bias)
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5 + bias;
            
            % Compute noisy range measurement
            epsR = self.measNoiseCovarSqrt * randn();
            
            R = R + epsR;
        end


        function R = computeMeasurement(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params(self.xIdx);
            ySensor = self.params(self.yIdx);
            zSensor = self.params(self.zIdx);
            bias = self.params(self.biasIdx);
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;
            
            % Compute range (with bias)
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5 + bias;
        end
    end
end
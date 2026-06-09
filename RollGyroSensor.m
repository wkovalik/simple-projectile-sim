classdef RollGyroSensor < Sensor
    properties (SetAccess = private)
        biasIdx = 0;
    end

    properties (Constant)
        nMeas = 1;

        DEFAULT_BIAS = 0;
    end


    methods
        % Constructor ==============================================================================

        function self = RollGyroSensor()
            self = self@Sensor();

            self.paramDefs.bias = ParamDef(self.DEFAULT_BIAS);

            self.updateParams();
        end


        % Update methods ===========================================================================

        function updateParams(self)
            self.params = [];

            self.biasIdx = 0;

            self.biasIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.bias.value];
        end


        function updateConsideredParams(self)
            self.consideredParams = [];
            self.consideredParamCovar = [];
            self.consideredParamIdxs = [];

            % bias
            if self.paramDefs.bias.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.bias.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.bias.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.biasIdx];
            end
        end


        % Measurement methods ======================================================================

        function p = takeMeasurement(self, state)
            % Get projectile state
            p = state(7);

            % Get sensor parameters
            bias = self.params(self.biasIdx);

            % Compute roll rate (with bias)
            p = p + bias;
            
            % Compute noisy roll rate measurement
            epsp = self.measNoiseCovarSqrt * randn();

            p = p + epsp;
        end


        function p = computeMeasurement(self, state)
            % Get projectile state
            p = state(7);

            % Get sensor parameters
            bias = self.params(self.biasIdx);

            % Compute roll rate (with bias)
            p = p + bias;
        end
    end
end
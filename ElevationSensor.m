classdef ElevationSensor < Sensor
    properties (Constant)
        N_MEASUREMENTS = 1;
    end

    methods
        % Constructor ==============================================================================
        function self = ElevationSensor()
            self = self@Sensor();

            self.params.x = Param();
            self.params.y = Param();
        end
        
        % Methods ==================================================================================
        function sampleMeasurement(self, time, state)
            x = state(1);
            y = state(2);
            
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            
            trueTheta = atan2(dySensor, dxSensor);
            
            noiseSigma = self.noiseCovar ^ 0.5;

            theta = trueTheta + noiseSigma * randn();
            
            self.sampleTime = time;
            self.sampledMeasurement = theta;
        end

        % Setters ==================================================================================
        % TODO
    end
end
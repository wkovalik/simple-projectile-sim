classdef RangeSensor < Sensor
    properties (Constant)
        N_MEASUREMENTS = 1;
    end

    methods
        % Constructor ==============================================================================
        function self = RangeSensor()
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
            
            trueR = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            noiseSigma = self.noiseCovar ^ 0.5;

            R = trueR + noiseSigma * randn();
            
            self.sampleTime = time;
            self.sampledMeasurement = R;
        end

        % Setters ==================================================================================
        % TODO
    end
end
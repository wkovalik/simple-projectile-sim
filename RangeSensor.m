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
            self.params.z = Param();
        end
        
        % Methods ==================================================================================
        function sampleMeasurement(self, time, state)
            x = state(1);
            y = state(2);
            z = state(3);
            
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            zSensor = self.params.z.value;
            
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;
            
            trueR = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
            
            noiseSigma = self.noiseCovar ^ 0.5;

            R = trueR + noiseSigma * randn();
            
            self.sampleTime = time;
            self.sampledMeasurement = R;
        end

        % Setters ==================================================================================
        % TODO
    end
end
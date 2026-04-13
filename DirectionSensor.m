classdef DirectionSensor < Sensor
    properties (Constant)
        N_MEASUREMENTS = 2;
    end

    methods
        % Constructor ==============================================================================
        function self = DirectionSensor()
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
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            trueAz = atan2(dySensor, dxSensor);
            trueEl = -atan2(dzSensor, drSensor);
            
            trueDir = [trueAz; trueEl];
            
            noiseSigma = chol(self.noiseCovar);

            dir = trueDir + noiseSigma * randn(self.N_MEASUREMENTS, 1);
            
            self.sampleTime = time;
            self.sampledMeasurement = dir;
        end

        % Setters ==================================================================================
        % TODO
    end
end
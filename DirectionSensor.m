classdef DirectionSensor < Sensor
    properties (Constant)
        N_MEASUREMENTS = 2;
    end

    methods
        % Constructor ==============================================================================
        function self = DirectionSensor()
            self = self@Sensor();

            self.params.x = Param();  % Position x-component (m)
            self.params.y = Param();  % Position y-component (m)
            self.params.z = Param();  % Position z-component (m)
        end
        

        % Methods ==================================================================================
        function sampleMeasurement(self, time, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            zSensor = self.params.z.value;
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;

            % Compute relative downrange
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            % Compute true azimuth and elevation
            trueAz = atan2(dySensor, dxSensor);
            trueEl = -atan2(dzSensor, drSensor);
            
            trueDir = [trueAz; trueEl];
            
            % Compute noisy azimuth and elevation measurements
            noiseSigma = chol(self.noiseCovar);  % Since covariance is a matrix here

            dir = trueDir + noiseSigma * randn(self.N_MEASUREMENTS, 1);
            
            % Store measurements
            self.sampleTime = time;
            self.sampledMeasurement = dir;
        end


        % Setters ==================================================================================
        % TODO
    end
end
classdef RangeSensor < Sensor
    properties (Constant)
        N_MEASUREMENTS = 1;
    end

    methods
        % Constructor ==============================================================================
        function self = RangeSensor()
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
            
            % Compute true range
            trueR = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
            
            % Compute noisy range measurement
            noiseSigma = self.noiseCovar ^ 0.5;
            
            R = trueR + noiseSigma * randn();
            
            % Store measurement
            self.sampleTime = time;
            self.sampledMeasurement = R;
        end


        % Setters ==================================================================================
        % TODO
    end
end
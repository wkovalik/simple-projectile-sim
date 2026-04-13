classdef Sensor < handle
    properties
        sensorID
        params

        samplePeriod
        sampleTime
        sampledMeasurement

        noiseCovar
    end

    properties (Dependent)
        invNoiseCovar

        nextSampleTime
    end
    
    properties (Abstract, Constant)
        N_MEASUREMENTS
    end

    methods
        % Constructor ==============================================================================
        function self = Sensor()
            self.samplePeriod = 0;   % (s)
            
            self.sampleTime = -999;  % (s) Not initialized to 0 to allow for a measurement at time = 0
            self.sampledMeasurement = zeros(self.N_MEASUREMENTS, 1);

            self.noiseCovar = zeros(self.N_MEASUREMENTS);
        end


        % Methods ==================================================================================
        function shouldSample = checkIfShouldSample(self, time)
            % Determines whether sensor is ready to take a new measurement
            
            shouldSample = (time >= (self.nextSampleTime - Constants.TIME_TOLERANCE));  % See Constants:Note 1
        end


        % Getters ==================================================================================
        function invNoiseCovar = get.invNoiseCovar(self)
            invNoiseCovar = inv(self.noiseCovar);
        end

        function nextSampleTime = get.nextSampleTime(self)
            nextSampleTime = self.sampleTime + self.samplePeriod;
        end

        
        % Setters ==================================================================================
        % TODO params
    end

    methods (Abstract)
        sampleMeasurement(self, time, state)
    end
end

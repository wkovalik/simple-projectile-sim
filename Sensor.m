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
    
    % TODO: Find out why this doesn't work
    % properties (Abstract)
    %     N_MEASUREMENTS
    % end

    methods
        % Constructor ==============================================================================
        function self = Sensor()
            self.samplePeriod = 0;
            
            self.sampleTime = -999;  % Not initialized to 0 to allow for a measurement at time = 0
            self.sampledMeasurement = 0;

            self.noiseCovar = 0;
        end

        % Methods ==================================================================================
        function shouldSample = checkIfShouldSample(self, time)
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
        % Methods ==================================================================================
        sampleMeasurement(self, time, state)
    end
end

classdef RollGyroSensor < Sensor
    properties
        projectileDynamics
    end

    properties (Constant)
        nMeas = 1;
    end

    methods
        % Constructor ==============================================================================

        function self = RollGyroSensor()
            self = self@Sensor();
        end


        % Measurement methods ======================================================================

        function p = takeMeasurement(self, state)
            % Get projectile state
            p = state(7);
            
            % Compute noisy acceleration measurement
            epsp = self.measNoiseCovarSqrt * randn();

            p = p + epsp;
        end


        function p = computeMeasurement(~, state)
            % Get projectile state
            p = state(7);
        end


        % Setters ==================================================================================

        function set.projectileDynamics(self, projectileDynamics)
            if Settings.VALIDATE_FLAG
                self.projectileDynamics = Validator.validateType(projectileDynamics, "ProjectileDynamics");
            else
                self.projectileDynamics = projectileDynamics;
            end
        end
    end
end
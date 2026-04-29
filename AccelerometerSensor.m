classdef AccelerometerSensor < Sensor
    properties
        projectileDynamics
    end

    properties (Constant)
        nMeas = 1;
    end

    methods
        % Constructor ==============================================================================

        function self = AccelerometerSensor(projectile, projectileDynamics)
            if nargin ~= 2
                error("Not enough input arguments. Requires projectile and projectileDynamics.")
            end
            
            self = self@Sensor(projectile, []);
            self.projectileDynamics = projectileDynamics;
        end


        % Measurement methods ======================================================================

        function a = takeMeasurement(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            % Get projectile parameters
            m = self.projectile.params(self.projectile.mIdx);
            
            % Compute acceleration
            FGrav = self.projectileDynamics.computeGravityForce();
            FAero = self.projectileDynamics.computeAeroForce(state);

            F = FGrav + FAero;
            a = F / m;
            
            % Compute component along velocity
            V = (vx ^ 2 + vy ^ 2 + vz ^ 2) ^ 0.5;
            unit_v = [vx; vy; vz] / V;

            a = dot(a, unit_v);
            
            % Compute noisy acceleration measurement
            epsa = self.measNoiseStdDev * randn();

            a = a + epsa;
        end


        function a = computeMeasurement(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            % Get projectile parameters
            m = self.projectile.params(self.projectile.mIdx);
            
            % Compute acceleration
            FGrav = self.projectileDynamics.computeGravityForce();
            FAero = self.projectileDynamics.computeAeroForce(state);

            F = FGrav + FAero;
            a = F / m;
            
            % Compute component along velocity
            V = (vx ^ 2 + vy ^ 2 + vz ^ 2) ^ 0.5;
            unit_v = [vx; vy; vz] / V;

            a = dot(a, unit_v);
        end
    end
end
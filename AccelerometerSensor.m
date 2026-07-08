classdef AccelerometerSensor < Sensor
    properties
        projectileDynamics
    end

    properties (SetAccess = private)
        biasIdx = 0;
    end

    properties (Constant)
        nMeas = 1;

        DEFAULT_BIAS = 0;
    end

    methods
        % Constructor ==============================================================================

        function self = AccelerometerSensor(projectile, projectileDynamics)
            if nargin ~= 2
                error("Not enough input arguments. Requires projectile and projectileDynamics.")
            end
            
            self = self@Sensor();

            self.projectile = projectile;
            self.projectileDynamics = projectileDynamics;

            self.paramDefs.bias = ParamDef(self.DEFAULT_BIAS);

            self.updateParams();
        end


        % Update methods ===========================================================================

        function updateParams(self)
            self.params = [];

            self.biasIdx = 0;

            self.biasIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.bias.value];
        end


        function updateConsideredParams(self)
            self.consideredParams = [];
            self.consideredParamCovar = [];
            self.consideredParamIdxs = [];

            % bias
            if self.paramDefs.bias.isConsidered
                self.consideredParams = [self.consideredParams; self.paramDefs.bias.value];
                self.consideredParamCovar = blkdiag(self.consideredParamCovar, self.paramDefs.bias.covar);
                self.consideredParamIdxs = [self.consideredParamIdxs; self.biasIdx];
            end
        end
        

        % Measurement methods ======================================================================

        function a = takeMeasurement(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            % Get projectile parameters
            m = self.projectile.params(self.projectile.mIdx);

            % Get sensor parameters
            bias = self.params(self.biasIdx);
            
            % Compute acceleration
            FGrav = self.projectileDynamics.computeGravityForce();
            FAero = self.projectileDynamics.computeAeroForce(state);

            F = FGrav + FAero;
            a = F / m;
            
            % Compute component along velocity (with bias)
            V = (vx ^ 2 + vy ^ 2 + vz ^ 2) ^ 0.5;
            unit_v = [vx; vy; vz] / V;

            a = dot(a, unit_v) + bias;
            
            % Compute noisy acceleration measurement
            epsa = self.measNoiseCovarSqrt * randn();

            a = a + epsa;
        end


        function a = computeMeasurement(self, state)
            % Get projectile state
            vx = state(4);
            vy = state(5);
            vz = state(6);
            
            % Get projectile parameters
            m = self.projectile.params(self.projectile.mIdx);

            % Get sensor parameters
            bias = self.params(self.biasIdx);
            
            % Compute acceleration
            FGrav = self.projectileDynamics.computeGravityForce();
            FAero = self.projectileDynamics.computeAeroForce(state);

            F = FGrav + FAero;
            a = F / m;
            
            % Compute component along velocity (with bias)
            V = (vx ^ 2 + vy ^ 2 + vz ^ 2) ^ 0.5;
            unit_v = [vx; vy; vz] / V;

            a = dot(a, unit_v) + bias;
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
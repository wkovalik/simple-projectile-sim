classdef DirectionSensor < Sensor
    properties (SetAccess = private)
        xIdx = 0;
        yIdx = 0;
        zIdx = 0;
    end

    properties (Constant)
        nMeas = 2;

        DEFAULT_X = 0;
        DEFAULT_Y = 0;
        DEFAULT_Z = 0;
    end

    methods
        % Constructor ==============================================================================

        function self = DirectionSensor(projectile, planet)
            if nargin == 0
                projectile = [];
                planet = [];
            end

            self = self@Sensor(projectile, planet);

            self.paramDefs.x = ParamDef(self.DEFAULT_X);
            self.paramDefs.y = ParamDef(self.DEFAULT_Y);
            self.paramDefs.z = ParamDef(self.DEFAULT_Z);

            self.updateParams();
        end

        
        % Update methods ===========================================================================

        function updateParams(self)
            self.params = [];

            self.xIdx = 0;
            self.yIdx = 0;
            self.zIdx = 0;

            self.xIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.x.value];

            self.yIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.y.value];

            self.zIdx = self.nParams + 1;
            self.params = [self.params; self.paramDefs.z.value];
        end
        

        % Measurement methods ======================================================================

        function dir = takeMeasurement(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params(self.xIdx);
            ySensor = self.params(self.yIdx);
            zSensor = self.params(self.zIdx);
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;

            % Compute true relative downrange
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            % Compute true azimuth and elevation
            az = atan2(dySensor, dxSensor);
            el = -atan2(dzSensor, drSensor);
            
            dir = [az; el];
            
            % Compute noisy azimuth and elevation measurements
            epsDir = self.measNoiseStdDev * randn(self.nMeas, 1);

            dir = dir + epsDir;
        end


        function dir = computeMeasurement(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params(self.xIdx);
            ySensor = self.params(self.yIdx);
            zSensor = self.params(self.zIdx);
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;

            % Compute relative downrange
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            % Compute azimuth and elevation
            az = atan2(dySensor, dxSensor);
            el = -atan2(dzSensor, drSensor);
            
            dir = [az; el];
        end
    end
end
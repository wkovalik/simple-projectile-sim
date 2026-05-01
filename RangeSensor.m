classdef RangeSensor < Sensor
    properties (SetAccess = private)
        xIdx = 0;
        yIdx = 0;
        zIdx = 0;
    end

    properties (Constant)
        nMeas = 1;

        DEFAULT_X = 0;
        DEFAULT_Y = 0;
        DEFAULT_Z = 0;
    end


    methods
        % Constructor ==============================================================================

        function self = RangeSensor()
            self = self@Sensor();

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
        
        function R = takeMeasurement(self, state)
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
            
            % Compute true range
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
            
            % Compute noisy range measurement
            epsR = self.measNoiseStdDev * randn();
            
            R = R + epsR;
        end


        function R = computeMeasurement(self, state)
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
            
            % Compute range
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
        end
    end
end
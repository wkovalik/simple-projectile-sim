classdef ElevationSensorModel < SensorModel
    properties (Constant)
        N_MEASUREMENTS = 1;
    end
    
    methods
        % Constructor ==============================================================================
        function self = ElevationSensorModel()
            self = self@SensorModel();

            self.params.x = Param();
            self.params.y = Param();
        end

        % Methods ==================================================================================
        function theta = computeMeasurement(self, state)
            x = state(1);
            y = state(2);
            
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            
            theta = atan2(dySensor, dxSensor);
        end

        function H = computeStateJacobian(self, state)
            x = state(1);
            y = state(2);
            
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            
            denom = 1 + (dySensor / dxSensor) ^ 2;
            
            H = zeros(self.N_MEASUREMENTS, self.N_PROJECTILE_STATES);
            
            H(1, 1) = -(dySensor / dxSensor ^ 2) / denom;
            H(1, 2) = (1 / dxSensor) / denom;
        end

        function H = computeDragJacobian(~, ~)
            H = 0;
        end

        function H = computeWindJacobian(~, ~)
            H = 0;
        end
    end
end

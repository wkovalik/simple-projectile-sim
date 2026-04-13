classdef RangeSensorModel < SensorModel
    properties (Constant)
        N_MEASUREMENTS = 1;
    end

    methods
        % Constructor ==============================================================================
        function self = RangeSensorModel()
            self = self@SensorModel();

            self.params.x = Param();
            self.params.y = Param();
            self.params.z = Param();
        end

        % Methods ==================================================================================
        function R = computeMeasurement(self, state)
            x = state(1);
            y = state(2);
            z = state(3);
            
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            zSensor = self.params.z.value;
            
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;
            
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
        end

        function H = computeStateJacobian(self, state)
            x = state(1);
            y = state(2);
            z = state(3);
            
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            zSensor = self.params.z.value;
            
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;
            
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
            
            H = zeros(self.N_MEASUREMENTS, self.N_PROJECTILE_STATES);
            
            H(1, 1) = dxSensor / R;
            H(1, 2) = dySensor / R;
            H(1, 3) = dzSensor / R;
        end

        function H = computeDragJacobian(~, ~)
            H = 0;
        end

        function H = computeWindJacobian(~, ~)
            H = 0;
        end
    end
end

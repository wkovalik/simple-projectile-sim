classdef DirectionSensorModel < SensorModel
    properties (Constant)
        N_MEASUREMENTS = 2;
    end
    
    methods
        % Constructor ==============================================================================
        function self = DirectionSensorModel()
            self = self@SensorModel();

            self.params.x = Param();
            self.params.y = Param();
            self.params.z = Param();
        end

        % Methods ==================================================================================
        function dir = computeMeasurement(self, state)
            x = state(1);
            y = state(2);
            z = state(3);
            
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            zSensor = self.params.z.value;
            
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;

            az = atan2(dySensor, dxSensor);
            el = -atan2(dzSensor, drSensor);

            dir = [az; el];
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
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            azDenom = 1 + (dySensor / dxSensor) ^ 2;
            elDenom = 1 + (dzSensor / drSensor) ^ 2;
            
            H = zeros(self.N_MEASUREMENTS, self.N_PROJECTILE_STATES);
            
            H(1, 1) = -(dySensor / dxSensor ^ 2) / azDenom;
            H(1, 2) = (1 / dxSensor) / azDenom;

            H(2, 1) = (dxSensor * dzSensor / drSensor ^ 3) / elDenom;
            H(2, 2) = (dySensor * dzSensor / drSensor ^ 3) / elDenom;
            H(2, 3) = -(1 / drSensor) / elDenom;
        end

        function H = computeDragJacobian(~, ~)
            H = [0; 0];
        end

        function H = computeWindJacobian(~, ~)
            H = [0; 0];
        end
    end
end

classdef DirectionSensorModel < SensorModel
    properties (Constant)
        N_MEASUREMENTS = 2;
    end
    
    methods
        % Constructor ==============================================================================
        function self = DirectionSensorModel()
            self = self@SensorModel();

            self.params.x = Param();  % Position x-component (m)
            self.params.y = Param();  % Position y-component (m)
            self.params.z = Param();  % Position z-component (m)
        end


        % Methods ==================================================================================
        function dir = computeMeasurement(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            zSensor = self.params.z.value;
            
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

        function H = computeStateJacobian(self, state)
            % Get projectile state
            x = state(1);
            y = state(2);
            z = state(3);
            
            % Get sensor parameters
            xSensor = self.params.x.value;
            ySensor = self.params.y.value;
            zSensor = self.params.z.value;
            
            % Compute relative position
            dxSensor = x - xSensor;
            dySensor = y - ySensor;
            dzSensor = z - zSensor;

            % Compute relative downrange
            drSensor = (dxSensor ^ 2 + dySensor ^ 2) ^ 0.5;
            
            % Compute denominators for partial derivatives
            azDenom = 1 + (dySensor / dxSensor) ^ 2;
            elDenom = 1 + (dzSensor / drSensor) ^ 2;
            
            % Build Jacobian -----------------------------------------------------------------------
            H = zeros(self.N_MEASUREMENTS, self.N_PROJECTILE_STATES);
            
            % Partial derivatives w.r.t. x
            H(1, 1) = -(dySensor / dxSensor ^ 2) / azDenom;
            H(2, 1) = (dxSensor * dzSensor / drSensor ^ 3) / elDenom;

            % Partial derivatives w.r.t. y
            H(1, 2) = (1 / dxSensor) / azDenom;
            H(2, 2) = (dySensor * dzSensor / drSensor ^ 3) / elDenom;

            % Partial derivatives w.r.t. z
            H(2, 3) = -(1 / drSensor) / elDenom;
        end

        function H = computeDragJacobian(self, ~)
            % Build Jacobian -----------------------------------------------------------------------

            % Partial derivatives w.r.t. CD
            H = zeros(self.N_MEASUREMENTS, 1);
        end

        function H = computeWindJacobian(self, ~)
            % Build Jacobian -----------------------------------------------------------------------

            % Partial derivatives w.r.t. vWindx
            H = zeros(self.N_MEASUREMENTS, 1);
        end
    end
end

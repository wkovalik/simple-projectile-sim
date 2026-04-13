classdef RangeSensorModel < SensorModel
    properties (Constant)
        N_MEASUREMENTS = 1;
    end

    methods
        % Constructor ==============================================================================
        function self = RangeSensorModel()
            self = self@SensorModel();

            self.params.x = Param();  % Position x-component (m)
            self.params.y = Param();  % Position y-component (m)
            self.params.z = Param();  % Position z-component (m)
        end


        % Methods ==================================================================================
        function R = computeMeasurement(self, state)
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
            
            % Compute range
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
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
            
            % Compute range
            R = (dxSensor ^ 2 + dySensor ^ 2 + dzSensor ^ 2) ^ 0.5;
            
            % Build Jacobian -----------------------------------------------------------------------
            H = zeros(self.N_MEASUREMENTS, self.N_PROJECTILE_STATES);
            
            % Partial derivatives w.r.t. x
            H(1, 1) = dxSensor / R;

            % Partial derivatives w.r.t. y
            H(1, 2) = dySensor / R;

            % Partial derivatives w.r.t. z
            H(1, 3) = dzSensor / R;
        end

        function H = computeDragJacobian(~, ~)
            % Build Jacobian -----------------------------------------------------------------------

            % Partial derivatives w.r.t. CD
            H = 0;
        end

        function H = computeWindxJacobian(~, ~)
            % Build Jacobian -----------------------------------------------------------------------
            
            % Partial derivatives w.r.t. vWindx
            H = 0;
        end

        function H = computeWindyJacobian(~, ~)
            % Build Jacobian -----------------------------------------------------------------------
            
            % Partial derivatives w.r.t. vWindy
            H = 0;
        end
    end
end

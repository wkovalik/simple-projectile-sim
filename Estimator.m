classdef Estimator < handle
    properties
        projectileModelDynamics

        projectileModel
        planetModel

        integrator
        propagator

        sensorModelArray
        sensorModelIDs
        sensorModelMap
    end

    properties (SetAccess = protected)
        includeParamSTM = false;
    end


    methods
        % Constructor ==============================================================================

        function self = Estimator(projectileModelDynamics, sensorModelArray, integrator)
            % Set handle for projectile model dynamics
            self.projectileModelDynamics = projectileModelDynamics;
            
            % Set handle for projectile and planet models (from dynamics model)
            self.projectileModel = self.projectileModelDynamics.projectile;
            self.planetModel = self.projectileModelDynamics.planet;
            
            % Set handle for integrator
            if nargin == 3
                self.integrator = integrator;
            elseif nargin == 2
                self.integrator = Integrator();
            else
                error("Not enough input arguments. Requires at least projectileModelDynamics and sensorModelArray.")
            end
            
            % Create propagator
            self.propagator = Propagator(self.projectileModelDynamics, self.integrator);
            
            % --------------------------------------------------------------------------------------
            
            % Store array of sensor models
            self.sensorModelArray = sensorModelArray;

            nSensorModels = length(self.sensorModelArray);
            
            % Build array of sensor IDs and dictionary of ID -> sensor model mappings
            self.sensorModelIDs = zeros(nSensorModels, 1);
            self.sensorModelMap = dictionary();

            for i = 1:nSensorModels
                % Pass handles for projectile and planet models to each sensor model (necessary to compute Jacobians)
                sensorModelArray{i}.projectile = self.projectileModel;
                sensorModelArray{i}.planet = self.planetModel;
                
                % Add each sensor model to ID array and dictionary
                self.sensorModelIDs(i) = sensorModelArray{i}.ID;
                self.sensorModelMap = insert(self.sensorModelMap, sensorModelArray{i}.ID, sensorModelArray(i));
            end
        end


        % Setters ==================================================================================

        function set.projectileModelDynamics(self, projectileModelDynamics)
            if Settings.VALIDATE_FLAG
                self.projectileModelDynamics = Validator.validateType(projectileModelDynamics, "ProjectileDynamics");
            else
                self.projectileModelDynamics = projectileModelDynamics;
            end
        end

        function set.projectileModel(self, projectileModel)
            if Settings.VALIDATE_FLAG
                self.projectileModel = Validator.validateType(projectileModel, "Projectile");
            else
                self.projectileModel = projectileModel;
            end
        end

        function set.planetModel(self, planetModel)
            if Settings.VALIDATE_FLAG
                self.planetModel = Validator.validateType(planetModel, "Planet");
            else
                self.planetModel = planetModel;
            end
        end

        function set.integrator(self, integrator)
            if Settings.VALIDATE_FLAG
                self.integrator = Validator.validateType(integrator, "Integrator");
            else
                self.integrator = integrator;
            end
        end

        function set.propagator(self, propagator)
            if Settings.VALIDATE_FLAG
                self.propagator = Validator.validateType(propagator, "Propagator");
            else
                self.propagator = propagator;
            end
        end

        function set.sensorModelArray(self, sensorModelArray)
            if Settings.VALIDATE_FLAG
                self.sensorModelArray = Validator.validateType(sensorModelArray, "cell");
            else
                self.sensorModelArray = sensorModelArray;
            end
        end

        function set.sensorModelIDs(self, sensorModelIDs)
            if Settings.VALIDATE_FLAG
                self.sensorModelIDs = Validator.validateType(sensorModelIDs, "double");
            else
                self.sensorModelIDs = sensorModelIDs;
            end
        end

        function set.sensorModelMap(self, sensorModelMap)
            if Settings.VALIDATE_FLAG
                self.sensorModelMap = Validator.validateType(sensorModelMap, "dictionary");
            else
                self.sensorModelMap = sensorModelMap;
            end
        end

        function set.includeParamSTM(self, includeParamSTM)
            if Settings.VALIDATE_FLAG
                self.includeParamSTM = Validator.validateType(includeParamSTM, "logical");
            else
                self.includeParamSTM = includeParamSTM;
            end
        end
    end

    methods (Abstract)
        % Solve method =============================================================================

        solve(self, measHistory);
    end
end
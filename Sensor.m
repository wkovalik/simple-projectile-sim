classdef Sensor < handle
    properties
        ID

        paramDefs

        params = [];

        % considerParams = [];
        % considerParamCovar = [];

        measNoiseCovar = [];
        invMeasNoiseCovar = [];
        measNoiseStdDev = [];

        samplePeriod = 0;
        nextSampleTime = 0;

        projectile
        planet
    end

    properties (SetAccess = private)
        nParams = 0;
        
        % nConsiderParams = 0;
        % considerParamIdxs = [];
    end
    
    properties (Abstract, Constant)
        nMeas
    end


    methods
        % Constructor ==============================================================================

        function self = Sensor(projectile, planet)
            if nargin == 2
                self.projectile = projectile;
                self.planet = planet;
            end

            self.measNoiseCovar = zeros(self.nMeas);
        end


        % Update methods ===========================================================================

        function updateParams(self)
            self.params = [];
        end

        
        % Jacobian methods =========================================================================

        function H = computeStateJacobian(self, state)
            nMeas = self.nMeas;
            nStates = self.projectile.nStates;

            H = zeros(nMeas, nStates);

            for i = 1:nStates
                delta = Settings.DEFAULT_JACOBIAN_PERT_FACTOR * (1 + abs(state(i)));

                statePlus = state;
                statePlus(i) = statePlus(i) + delta;

                measurementPlus = self.computeMeasurement(statePlus);

                stateMinus = state;
                stateMinus(i) = stateMinus(i) - delta;

                measurementMinus = self.computeMeasurement(stateMinus);

                H(:, i) = (measurementPlus - measurementMinus) / (2 * delta);
            end
        end


        function H = computeParamJacobian(self, state)
            nMeas = self.nMeas;
            nEstimatedProjectileParams = self.projectile.nEstimatedParams;
            nEstimatedPlanetParams = self.planet.nEstimatedParams;
            nEstimatedParams = nEstimatedProjectileParams + nEstimatedPlanetParams;

            H = zeros(nMeas, nEstimatedParams);
            
            for i = 1:nEstimatedProjectileParams
                paramIdx = self.projectile.estimatedParamIdxs(i);
                param = self.projectile.params(paramIdx);

                delta = Settings.DEFAULT_JACOBIAN_PERT_FACTOR * (1 + abs(param));

                paramPlus = param + delta;
                self.projectile.params(paramIdx) = paramPlus;

                measurementPlus = self.computeMeasurement(state);

                paramMinus = param - delta;
                self.projectile.params(paramIdx) = paramMinus;

                measurementMinus = self.computeMeasurement(state);

                H(:, i) = (measurementPlus - measurementMinus) / (2 * delta);

                self.projectile.params(paramIdx) = param;
            end

            for i = 1:nEstimatedPlanetParams
                paramIdx = self.planet.estimatedParamIdxs(i);
                param = self.planet.params(paramIdx);

                delta = Settings.DEFAULT_JACOBIAN_PERT_FACTOR * (1 + abs(param));

                paramPlus = param + delta;
                self.planet.params(paramIdx) = paramPlus;

                measurementPlus = self.computeMeasurement(state);

                paramMinus = param - delta;
                self.planet.params(paramIdx) = paramMinus;

                measurementMinus = self.computeMeasurement(state);

                H(:, nEstimatedProjectileParams + i) = (measurementPlus - measurementMinus) / (2 * delta);

                self.planet.params(paramIdx) = param;
            end
        end


        % Helper methods ===========================================================================

        function isReadyToSample = shouldTakeMeasurement(self, time)
            % Determines whether sensor is ready to take a new measurement
            
            isReadyToSample = (time >= (self.nextSampleTime - Settings.DEFAULT_TIME_TOL));  % See Settings:Note 1

            if isReadyToSample
                self.nextSampleTime = time + self.samplePeriod;
            end
        end

        
        % Setters ==================================================================================

        function set.ID(self, ID)
            if Settings.VALIDATE_FLAG
                self.ID = Validator.validateType(ID, "double");
            else
                self.ID = ID;
            end
        end

        function set.paramDefs(self, paramDefs)
            if Settings.VALIDATE_FLAG
                self.paramDefs = Validator.validateFieldTypes(paramDefs, "ParamDef");
            else
                self.paramDefs = paramDefs;
            end
        end

        function set.params(self, params)
            if Settings.VALIDATE_FLAG
                self.params = Validator.validateType(params, "double");
            else
                self.params = params;
            end

            self.nParams = length(params);
        end

        function set.measNoiseCovar(self, measNoiseCovar)
            if Settings.VALIDATE_FLAG
                measNoiseCovar = Validator.validateType(measNoiseCovar, "double");
                self.measNoiseCovar = Validator.validateSize(measNoiseCovar, [self.nMeas, self.nMeas]);
            else
                self.measNoiseCovar = measNoiseCovar;
            end
            
            % Covariance must be nonsingular for these to be defined
            if det(measNoiseCovar) ~= 0
                self.invMeasNoiseCovar = inv(measNoiseCovar);
                self.measNoiseStdDev = chol(measNoiseCovar);   % Matrix square root
            end
        end

        function set.samplePeriod(self, samplePeriod)
            if Settings.VALIDATE_FLAG
                self.samplePeriod = Validator.validateType(samplePeriod, "double");
            else
                self.samplePeriod = samplePeriod;
            end
        end

        function set.nextSampleTime(self, nextSampleTime)
            if Settings.VALIDATE_FLAG
                self.nextSampleTime = Validator.validateType(nextSampleTime, "double");
            else
                self.nextSampleTime = nextSampleTime;
            end
        end
    end

    methods (Abstract)
        % Measurement methods ======================================================================

        takeMeasurement(self, state)
        computeMeasurement(self, state)
    end
end

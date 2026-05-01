classdef BatchEstimator < handle
    % TODO: Currently assumes init projectile time = estimate time epoch. Need to propagate if not

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

    properties (SetAccess = private)
        includeParamSTM = false;
    end


    methods
        % Constructor ==============================================================================

        function self = BatchEstimator(projectileModelDynamics, sensorModelArray, integrator)
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


        % Solve method =============================================================================

        function output = solve(self, measHistory)
            % Remove any measurements from sensors not included in sensor model array
            measHistory = measHistory(:, ismember(measHistory(1, :), self.sensorModelIDs));

            nSamples = size(measHistory, 2);
            
            measTimeHistory = measHistory(2, :);   % Get measurement time history
            finalMeasTime = measTimeHistory(end);  % Get time of final measurement

            % --------------------------------------------------------------------------------------
            nStates = self.projectileModel.nStates;
            nEstimatedProjectileParams = self.projectileModel.nEstimatedParams;
            nEstimatedPlanetParams = self.planetModel.nEstimatedParams;
            nEstimatedParams = nEstimatedProjectileParams + nEstimatedPlanetParams;

            self.includeParamSTM = logical(nEstimatedParams);
            self.projectileModelDynamics.includeParamSTM = logical(nEstimatedParams);

            nAugStates = nStates + nEstimatedParams;

            % Build prefit vectors and covariance matrices
            priorState = self.projectileModel.state;
            priorStateCovar = self.projectileModel.stateCovar;

            priorParams = [self.projectileModel.estimatedParams; self.planetModel.estimatedParams];
            priorParamCovar = blkdiag(self.projectileModel.estimatedParamCovar, self.planetModel.estimatedParamCovar);

            priorAugState = [priorState; priorParams];
            priorAugStateCovar = blkdiag(priorStateCovar, priorParamCovar);
            priorAugStateDelta = zeros(nStates + nEstimatedParams, 1);

            fprintf("Iteration\tEstimated State\n")
            fprintf("%i\t\t\t", 0)
            fprintf("%.4f\t", priorAugState(:))
            fprintf("\n")
            
            % Initialize output data structure
            nMaxIterations = Settings.DEFAULT_MAX_ITERS;

            output.iterationData = cell(1, nMaxIterations + 1);

            output.stateHistory = zeros(nStates, nMaxIterations + 1);
            output.stateCovarHistory = zeros(nStates ^ 2, nMaxIterations + 1);

            output.paramHistory = zeros(nEstimatedParams, nMaxIterations + 1);
            output.paramCovarHistory = zeros(nEstimatedParams ^ 2, nMaxIterations + 1);

            output.augStateHistory = zeros(nAugStates, nMaxIterations + 1);
            output.augStateCovarHistory = zeros(nAugStates ^ 2, nMaxIterations + 1);
            
            % Add prefit vectors and covariances to output
            output.stateHistory(:, 1) = priorState;
            output.stateCovarHistory(:, 1) = priorStateCovar(:);

            output.paramHistory(:, 1) = priorParams;
            output.paramCovarHistory(:, 1) = priorParamCovar(:);

            output.augStateHistory(:, 1) = priorAugState;
            output.augStateCovarHistory(:, 1) = priorAugStateCovar(:);

            hasConverged = false;

            % --------------------------------------------------------------------------------------
            % Begin estimation loop
            % --------------------------------------------------------------------------------------

            for ii = 1:nMaxIterations
                % Propagate nominal trajectory and STMs
                [nomTimeHistory, nomStateHistory, nomStateSTMHistory, nomParamSTMHistory] = ...
                    self.propagator.propagateWithSTM(finalMeasTime);
                
                output.iterationData{ii}.nomTimeHistory = nomTimeHistory;
                output.iterationData{ii}.nomStateHistory = nomStateHistory;
                
                % Resample nominal trajectory at measurement times (should be exact)
                nomStateHistory = Utils.resampleStateHistory(nomTimeHistory, nomStateHistory, measTimeHistory);
                nomStateSTMHistory = Utils.resampleStateHistory(nomTimeHistory, nomStateSTMHistory, measTimeHistory);
                if self.includeParamSTM
                    nomParamSTMHistory = Utils.resampleStateHistory(nomTimeHistory, nomParamSTMHistory, measTimeHistory);
                end
                
                % ----------------------------------------------------------------------------------
                
                % Initialize postfit normal vector and information matrix (i.e., inverse covariance matrix)
                postAugNormal = priorAugStateCovar \ priorAugStateDelta;
                postAugStateInvCovar = inv(priorAugStateCovar);
                
                % Initialize measurement residual history
                measResidualHistory = zeros(size(measHistory, 1), nSamples);
                measResidualHistory(1, :) = measHistory(1, :);
                measResidualHistory(2, :) = measHistory(2, :);
    
                for i = 1:nSamples
                    % Get sensor model for current measurement
                    sensorID = measHistory(1, i);
                    sensorModel = self.sensorModelMap{sensorID};
                    
                    % Get nominal state and STMs at current measurement time
                    nomState = nomStateHistory(:, i);
                    nomStateSTM = nomStateSTMHistory(:, i);
                    nomStateSTM = reshape(nomStateSTM, [nStates, nStates]);
                    if self.includeParamSTM
                        nomParamSTM = nomParamSTMHistory(:, i);
                        nomParamSTM = reshape(nomParamSTM, [nStates, nEstimatedParams]);
                    end
                    
                    % Get observed measurement and computed measurement at current measurement time
                    nMeas = sensorModel.nMeas;
                    iMeasEnd = 3 + (nMeas - 1);

                    observedMeas = measHistory(3:iMeasEnd, i);
                    computedMeas = sensorModel.computeMeasurement(nomState);
    
                    % Compute measurement residual
                    measResidual = observedMeas - computedMeas;
                    measResidualHistory(3:iMeasEnd, i) = measResidual;
                    
                    % Compute measurement sensitivity matrices (i.e., Jacobians) and map to initial time
                    stateH = sensorModel.computeStateJacobian(nomState);
                    mappedStateH = stateH * nomStateSTM;

                    if self.includeParamSTM
                        paramH = sensorModel.computeParamJacobian(nomState);
                        mappedParamH = stateH * nomParamSTM + paramH;
                    else
                        mappedParamH = [];
                    end
                    
                    mappedH = [mappedStateH, mappedParamH];
                    
                    % Accumulate postfit normal vector and information matrix
                    invMeasNoiseCovar = sensorModel.invMeasNoiseCovar;
                    
                    nextAugNormal = mappedH' * invMeasNoiseCovar * measResidual;
                    nextAugStateInvCovar = mappedH' * invMeasNoiseCovar * mappedH;
    
                    postAugNormal = postAugNormal + nextAugNormal;
                    postAugStateInvCovar = postAugStateInvCovar + nextAugStateInvCovar;
                end

                output.iterationData{ii}.measResidualHistory = measResidualHistory;
                
                % ----------------------------------------------------------------------------------
                
                % Compute postfit state deviation and covariance
                postAugStateDelta = postAugStateInvCovar \ postAugNormal;
                postAugStateCovar = inv(postAugStateInvCovar);

                if ii == 1
                    % Determine if state has converged
                    if max(abs(postAugStateDelta ./ priorAugState)) < Settings.DEFAULT_CONVERGENCE_TOL
                        hasConverged = true;
                    end
                    
                    % Update prefit state (now becomes postfit state)
                    postAugState = priorAugState + postAugStateDelta;
                else
                    % Determine if state has converged
                    if max(abs(postAugStateDelta ./ postAugState)) < Settings.DEFAULT_CONVERGENCE_TOL
                        hasConverged = true;
                    end
                    
                    % Update postfit state
                    postAugState = postAugState + postAugStateDelta;
                end

                fprintf("%i\t\t\t", ii)
                fprintf("%.4f\t", postAugState(:))
                fprintf("\n")
                
                output.augStateHistory(:, ii + 1) = postAugState;
                output.augStateCovarHistory(:, ii + 1) = postAugStateCovar(:);
                
                % Shift prefit state deviation
                priorAugStateDelta = priorAugStateDelta - postAugStateDelta;

                % ----------------------------------------------------------------------------------
                
                % Extract postfit projectile state
                postState = postAugState(1:nStates);
                postStateCovar = postAugStateCovar(1:nStates, 1:nStates);

                output.stateHistory(:, ii + 1) = postState;
                output.stateCovarHistory(:, ii + 1) = postStateCovar(:);
                
                % Update projectile model state
                self.projectileModel.time = 0;  % See TODO
                self.projectileModel.state = postState;

                if self.includeParamSTM
                    % Extract postfit parameters
                    postParams = postAugState((nStates + 1):end);
                    postParamCovar = postAugStateCovar((nStates + 1):end, (nStates + 1):end);
    
                    output.paramHistory(:, ii + 1) = postParams;
                    output.paramCovarHistory(:, ii + 1) = postParamCovar(:);
    
                    % Extract postfit projectile parameters
                    postProjectileParams = postParams(1:nEstimatedProjectileParams);
                    postProjectileParamCovar = postParamCovar(1:nEstimatedProjectileParams, 1:nEstimatedProjectileParams);
                    
                    % Update projectile model parameters
                    self.projectileModel.params(self.projectileModel.estimatedParamIdxs) = postProjectileParams;
                    self.projectileModel.estimatedParams = postProjectileParams;
                    self.projectileModel.estimatedParamCovar = postProjectileParamCovar;
    
                    % Extract postfit planet parameters
                    postPlanetParams = postParams((nEstimatedProjectileParams + 1):end);
                    postPlanetParamCovar = postParamCovar((nEstimatedProjectileParams + 1):end, (nEstimatedProjectileParams + 1):end);
                    
                    % Update planet model parameters
                    self.planetModel.params(self.planetModel.estimatedParamIdxs) = postPlanetParams;
                    self.planetModel.estimatedParams = postPlanetParams;
                    self.planetModel.estimatedParamCovar = postPlanetParamCovar;
                end

                % ----------------------------------------------------------------------------------

                if hasConverged
                    fprintf("Converged!\n")  % Break out if converged
                    break

                elseif ii == nMaxIterations
                    warning("Failed to converge within maximum number of iterations.")
                end
            end
            
            % --------------------------------------------------------------------------------------
            % End estimation loop
            % --------------------------------------------------------------------------------------

            nIterations = ii;
            
            % Propagate postfit nominal trajectory
            [nomTimeHistory, nomStateHistory] = self.propagator.propagate(finalMeasTime);
                
            output.iterationData{nIterations + 1}.nomTimeHistory = nomTimeHistory;
            output.iterationData{nIterations + 1}.nomStateHistory = nomStateHistory;
            
            % Resample postfit nominal trajectory at measurement times (should be exact)
            nomStateHistory = Utils.resampleStateHistory(nomTimeHistory, nomStateHistory, measTimeHistory);
            
            % Initialize postfit measurement residual history
            measResidualHistory = zeros(size(measHistory, 1), nSamples);
            measResidualHistory(1, :) = measHistory(1, :);
            measResidualHistory(2, :) = measHistory(2, :);

            for i = 1:nSamples
                % Get sensor model for current measurement
                sensorID = measHistory(1, i);
                sensorModel = self.sensorModelMap{sensorID};
                
                % Get nominal state at current measurement time
                nomState = nomStateHistory(:, i);
                
                % Get observed measurement and computed measurement at current measurement time
                nMeas = sensorModel.nMeas;
                iMeasEnd = 3 + (nMeas - 1);

                observedMeas = measHistory(3:iMeasEnd, i);
                computedMeas = sensorModel.computeMeasurement(nomState);
                
                % Compute postfit measurement residual
                measResidual = observedMeas - computedMeas;
                measResidualHistory(3:iMeasEnd, i) = measResidual;
            end

            output.iterationData{nIterations + 1}.measResidualHistory = measResidualHistory;
            
            % Remove all unused entries in output data if converged early
            if hasConverged
                output.iterationData((nIterations + 2):end) = [];

                output.stateHistory(:, (nIterations + 2):end) = [];
                output.stateCovarHistory(:, (nIterations + 2):end) = [];

                output.paramHistory(:, (nIterations + 2):end) = [];
                output.paramCovarHistory(:, (nIterations + 2):end) = [];

                output.augStateHistory(:, (nIterations + 2):end) = [];
                output.augStateCovarHistory(:, (nIterations + 2):end) = [];
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
end

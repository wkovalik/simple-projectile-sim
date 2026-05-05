classdef BatchEstimator < Estimator
    % TODO: Currently assumes init projectile time = estimate time epoch. Need to propagate if not

    methods
        % Constructor ==============================================================================

        function self = BatchEstimator(varargin)
            self = self@Estimator(varargin{:});
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

            % Build prefit vectors and covariances at initial time
            priorState_0 = self.projectileModel.state;
            priorStateCovar_0 = self.projectileModel.stateCovar;

            priorParams = [self.projectileModel.estimatedParams; self.planetModel.estimatedParams];
            priorParamCovar = blkdiag(self.projectileModel.estimatedParamCovar, self.planetModel.estimatedParamCovar);

            priorAugState_0 = [priorState_0; priorParams];
            priorAugStateCovar_0 = blkdiag(priorStateCovar_0, priorParamCovar);
            priorAugStateDelta_0 = zeros(nAugStates, 1);

            fprintf("Iteration\tEstimated State\n")
            fprintf("%i\t\t\t", 0)
            fprintf("%.4f\t", priorAugState_0(:))
            fprintf("\n")
            
            % Initialize output data structure
            nMaxIterations = Settings.DEFAULT_MAX_ITERS;

            output.perIterationData = cell(1, nMaxIterations + 1);

            output.iterations.state0 = zeros(nStates, nMaxIterations + 1);
            output.iterations.stateCovar0 = zeros(nStates ^ 2, nMaxIterations + 1);

            output.iterations.params = zeros(nEstimatedParams, nMaxIterations + 1);
            output.iterations.paramCovar = zeros(nEstimatedParams ^ 2, nMaxIterations + 1);

            output.iterations.augState0 = zeros(nAugStates, nMaxIterations + 1);
            output.iterations.augStateCovar0 = zeros(nAugStates ^ 2, nMaxIterations + 1);
            
            % Add prefit vectors and covariances to output
            output.iterations.state0(:, 1) = priorState_0;
            output.iterations.stateCovar0(:, 1) = priorStateCovar_0(:);

            output.iterations.params(:, 1) = priorParams;
            output.iterations.paramCovar(:, 1) = priorParamCovar(:);

            output.iterations.augState0(:, 1) = priorAugState_0;
            output.iterations.augStateCovar0(:, 1) = priorAugStateCovar_0(:);

            hasConverged = false;

            % --------------------------------------------------------------------------------------
            % Begin estimation loop
            % --------------------------------------------------------------------------------------

            for ii = 1:nMaxIterations
                % Propagate nominal trajectory and STMs
                [nomTimeHistory, nomStateHistory, stateSTMHistory, paramSTMHistory] = ...
                    self.propagator.propagateWithSTM(finalMeasTime);
                
                output.perIterationData{ii}.nomTimeHistory = nomTimeHistory;
                output.perIterationData{ii}.nomStateHistory = nomStateHistory;
                
                % Resample nominal trajectory at measurement times (should be exact)
                nomStateHistory = Utils.resampleStateHistory(nomTimeHistory, nomStateHistory, measTimeHistory);
                stateSTMHistory = Utils.resampleStateHistory(nomTimeHistory, stateSTMHistory, measTimeHistory);
                if self.includeParamSTM
                    paramSTMHistory = Utils.resampleStateHistory(nomTimeHistory, paramSTMHistory, measTimeHistory);
                end
                
                % ----------------------------------------------------------------------------------
                
                % Initialize postfit normal vector and information matrix (i.e., inverse covariance) for initial time
                postAugNormal_0 = priorAugStateCovar_0 \ priorAugStateDelta_0;
                postAugStateInvCovar_0 = inv(priorAugStateCovar_0);
                
                % Initialize measurement residual history
                measResidualHistory = zeros(size(measHistory, 1), nSamples);
                measResidualHistory(1, :) = measHistory(1, :);
                measResidualHistory(2, :) = measHistory(2, :);
    
                for i = 1:nSamples
                    % Get sensor model for current measurement
                    sensorID_i = measHistory(1, i);
                    sensorModel_i = self.sensorModelMap{sensorID_i};
                    
                    % Get nominal state and STMs at current measurement time
                    nomState_i = nomStateHistory(:, i);

                    stateSTM_i0 = stateSTMHistory(:, i);
                    stateSTM_i0 = reshape(stateSTM_i0, [nStates, nStates]);
                    if self.includeParamSTM
                        paramSTM_i0 = paramSTMHistory(:, i);
                        paramSTM_i0 = reshape(paramSTM_i0, [nStates, nEstimatedParams]);
                    end
                    
                    % Get observed measurement and computed measurement at current measurement time
                    nMeas = sensorModel_i.nMeas;
                    iMeasEnd = 3 + (nMeas - 1);

                    observedMeas_i = measHistory(3:iMeasEnd, i);
                    computedMeas_i = sensorModel_i.computeMeasurement(nomState_i);
    
                    % Compute measurement residual
                    measResidual_i = observedMeas_i - computedMeas_i;
                    measResidualHistory(3:iMeasEnd, i) = measResidual_i;
                    
                    % Compute measurement sensitivity matrices (i.e., Jacobians) at current measurement time
                    stateH_i = sensorModel_i.computeStateJacobian(nomState_i);
                    if self.includeParamSTM
                        paramH_i = sensorModel_i.computeParamJacobian(nomState_i);
                    end
                    
                    % Map measurement sensitivity matrices to initial time
                    mappedStateH_i0 = stateH_i * stateSTM_i0;
                    if self.includeParamSTM
                        mappedParamH_i0 = stateH_i * paramSTM_i0 + paramH_i;

                        mappedH_i0 = [mappedStateH_i0, mappedParamH_i0];
                    else
                        mappedH_i0 = mappedStateH_i0;
                    end
                    
                    % Get measurement noise covariance inverse
                    invMeasNoiseCovar_i = sensorModel_i.invMeasNoiseCovar;
                    
                    % Accumulate postfit normal vector and information matrix
                    addAugNormal_0 = mappedH_i0' * invMeasNoiseCovar_i * measResidual_i;
                    addAugStateInvCovar_0 = mappedH_i0' * invMeasNoiseCovar_i * mappedH_i0;
    
                    postAugNormal_0 = postAugNormal_0 + addAugNormal_0;
                    postAugStateInvCovar_0 = postAugStateInvCovar_0 + addAugStateInvCovar_0;
                end

                output.perIterationData{ii}.measResidualHistory = measResidualHistory;
                
                % ----------------------------------------------------------------------------------
                
                % Compute postfit state deviation and covariance at initial time
                postAugStateDelta_0 = postAugStateInvCovar_0 \ postAugNormal_0;
                postAugStateCovar_0 = inv(postAugStateInvCovar_0);

                if ii == 1
                    % Determine if state has converged (using norm convergence: ||x_i - x_(i-1)|| / ||x_(i-1)||)
                    if (norm(postAugStateDelta_0) / norm(priorAugState_0)) < Settings.DEFAULT_CONVERGENCE_TOL
                        hasConverged = true;
                    end
                    
                    % Update prefit state (now becomes postfit state)
                    postAugState_0 = priorAugState_0 + postAugStateDelta_0;
                else
                    % Determine if state has converged (using norm convergence: ||x_i - x_(i-1)|| / ||x_(i-1)||)
                    if (norm(postAugStateDelta_0) / norm(postAugState_0)) < Settings.DEFAULT_CONVERGENCE_TOL
                        hasConverged = true;
                    end
                    
                    % Update postfit state
                    postAugState_0 = postAugState_0 + postAugStateDelta_0;
                end

                fprintf("%i\t\t\t", ii)
                fprintf("%.4f\t", postAugState_0(:))
                fprintf("\n")
                
                output.iterations.augState0(:, ii + 1) = postAugState_0;
                output.iterations.augStateCovar0(:, ii + 1) = postAugStateCovar_0(:);
                
                % Shift prefit state deviation at initial time
                priorAugStateDelta_0 = priorAugStateDelta_0 - postAugStateDelta_0;

                % ----------------------------------------------------------------------------------
                
                % Extract postfit projectile state at initial time
                postState_0 = postAugState_0(1:nStates);
                postStateCovar_0 = postAugStateCovar_0(1:nStates, 1:nStates);

                output.iterations.state0(:, ii + 1) = postState_0;
                output.iterations.stateCovar0(:, ii + 1) = postStateCovar_0(:);
                
                % Update projectile model state (for nominal trajectory on next iteration)
                self.projectileModel.time = 0;  % See TODO
                self.projectileModel.state = postState_0;

                if self.includeParamSTM
                    % Extract postfit parameters
                    postParams = postAugState_0((nStates + 1):end);
                    postParamCovar = postAugStateCovar_0((nStates + 1):end, (nStates + 1):end);
    
                    output.iterations.params(:, ii + 1) = postParams;
                    output.iterations.paramCovar(:, ii + 1) = postParamCovar(:);
    
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
                
            output.perIterationData{nIterations + 1}.nomTimeHistory = nomTimeHistory;
            output.perIterationData{nIterations + 1}.nomStateHistory = nomStateHistory;
            
            % Resample postfit nominal trajectory at measurement times (should be exact)
            nomStateHistory = Utils.resampleStateHistory(nomTimeHistory, nomStateHistory, measTimeHistory);
            
            % Initialize postfit measurement residual history
            measResidualHistory = zeros(size(measHistory, 1), nSamples);
            measResidualHistory(1, :) = measHistory(1, :);
            measResidualHistory(2, :) = measHistory(2, :);

            for i = 1:nSamples
                % Get sensor model for current measurement
                sensorID_i = measHistory(1, i);
                sensorModel_i = self.sensorModelMap{sensorID_i};
                
                % Get nominal state at current measurement time
                nomState_i = nomStateHistory(:, i);
                
                % Get observed measurement and computed measurement at current measurement time
                nMeas = sensorModel_i.nMeas;
                iMeasEnd = 3 + (nMeas - 1);

                observedMeas_i = measHistory(3:iMeasEnd, i);
                computedMeas_i = sensorModel_i.computeMeasurement(nomState_i);
                
                % Compute postfit measurement residual
                measResidual_i = observedMeas_i - computedMeas_i;
                measResidualHistory(3:iMeasEnd, i) = measResidual_i;
            end

            output.perIterationData{nIterations + 1}.measResidualHistory = measResidualHistory;
            
            % Remove all unused entries in output data if converged early
            if hasConverged
                output.perIterationData((nIterations + 2):end) = [];

                output.iterations.state0(:, (nIterations + 2):end) = [];
                output.iterations.stateCovar0(:, (nIterations + 2):end) = [];

                output.iterations.params(:, (nIterations + 2):end) = [];
                output.iterations.paramCovar(:, (nIterations + 2):end) = [];

                output.iterations.augState0(:, (nIterations + 2):end) = [];
                output.iterations.augStateCovar0(:, (nIterations + 2):end) = [];
            end
        end
    end
end

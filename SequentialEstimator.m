classdef SequentialEstimator < Estimator
    % TODO: Currently assumes init projectile time = estimate time epoch. Need to propagate if not

    methods
        % Constructor ==============================================================================

        function self = SequentialEstimator(varargin)
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

            output.iterationData = cell(1, nMaxIterations + 1);

            output.state_0Iterations = zeros(nStates, nMaxIterations + 1);
            output.stateCovar_0Iterations = zeros(nStates ^ 2, nMaxIterations + 1);

            output.paramIterations = zeros(nEstimatedParams, nMaxIterations + 1);
            output.paramCovarIterations = zeros(nEstimatedParams ^ 2, nMaxIterations + 1);

            output.augState_0Iterations = zeros(nAugStates, nMaxIterations + 1);
            output.augStateCovar_0Iterations = zeros(nAugStates ^ 2, nMaxIterations + 1);
            
            % Add prefit vectors and covariances to output
            output.state_0Iterations(:, 1) = priorState_0;
            output.stateCovar_0Iterations(:, 1) = priorStateCovar_0(:);

            output.paramIterations(:, 1) = priorParams;
            output.paramCovarIterations(:, 1) = priorParamCovar(:);

            output.augState_0Iterations(:, 1) = priorAugState_0;
            output.augStateCovar_0Iterations(:, 1) = priorAugStateCovar_0(:);

            hasConverged = false;

            % --------------------------------------------------------------------------------------
            % Begin estimation loop
            % --------------------------------------------------------------------------------------

            for ii = 1:nMaxIterations
                % Propagate nominal trajectory and STMs
                [nomTimeHistory, nomStateHistory, stateSTMHistory, paramSTMHistory] = ...
                    self.propagator.propagateWithSTM(finalMeasTime);
                
                output.iterationData{ii}.nomTimeHistory = nomTimeHistory;
                output.iterationData{ii}.nomStateHistory = nomStateHistory;
                
                % Resample nominal trajectory at measurement times (should be exact)
                nomStateHistory = Utils.resampleStateHistory(nomTimeHistory, nomStateHistory, measTimeHistory);
                stateSTMHistory = Utils.resampleStateHistory(nomTimeHistory, stateSTMHistory, measTimeHistory);
                if self.includeParamSTM
                    paramSTMHistory = Utils.resampleStateHistory(nomTimeHistory, paramSTMHistory, measTimeHistory);
                end
                
                % ----------------------------------------------------------------------------------
                
                % Note: t_i = current measurement time
                %       t_j = t_(i-1) = previous measurement time
                
                % Initialize postfit state deviation and covariance at previous measurement time (i.e., initial time here)
                postAugStateDelta_j = priorAugStateDelta_0;
                postAugStateCovar_j = priorAugStateCovar_0;
                
                % Initialize STMs at previous measurement time (i.e., initial time here)
                invStateSTM_j0 = eye(nStates);
                if self.includeParamSTM
                    paramSTM_j0 = zeros(nStates, nEstimatedParams);
    
                    invSTM_j0 = [invStateSTM_j0,                   -invStateSTM_j0 * paramSTM_j0;
                                 zeros(nEstimatedParams, nStates),  eye(nEstimatedParams)];
                else
                    invSTM_j0 = invStateSTM_j0;
                end
    
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
    
                    if self.includeParamSTM
                        STM_i0 = [stateSTM_i0,                      paramSTM_i0;
                                  zeros(nEstimatedParams, nStates), eye(nEstimatedParams)];
                    else
                        STM_i0 = stateSTM_i0;
                    end
                    
                    % Compute step STM from previous measurement time to current measurement time
                    STM_ij = STM_i0 * invSTM_j0;
    
                    % Compute prefit state deviation and covariance at current measurement time
                    priorAugStateDelta_i = STM_ij * postAugStateDelta_j;
                    priorAugStateCovar_i = STM_ij * postAugStateCovar_j * STM_ij';
                    
                    % Above: Time update
                    % ------------------------------------------------------------------------------
                    % Below: Measurement update

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
    
                        H_i = [stateH_i, paramH_i];
                    else
                        H_i = stateH_i;
                    end
    
                    measNoiseCovar_i = sensorModel_i.measNoiseCovar;

                    % Compute filter gain (i.e., Kalman gain) matrix
                    measResidualGain_i = priorAugStateCovar_i * H_i' / (H_i * priorAugStateCovar_i * H_i' + measNoiseCovar_i);
    
                    % Compute postfit state deviation and covariance at current measurement time
                    postAugStateDelta_i = priorAugStateDelta_i + measResidualGain_i * (measResidual_i - H_i * priorAugStateDelta_i);
                    postAugStateCovar_i = (eye(nAugStates) - measResidualGain_i * H_i) * priorAugStateCovar_i * (eye(nAugStates) - measResidualGain_i * H_i)' + ...
                                          measResidualGain_i * measNoiseCovar_i * measResidualGain_i';

                    % Store results (current measurement time now becomes previous measurement time)
                    postAugStateDelta_j = postAugStateDelta_i;
                    postAugStateCovar_j = postAugStateCovar_i;
                    
                    invStateSTM_j0 = inv(stateSTM_i0);
                    if self.includeParamSTM
                        paramSTM_j0 = paramSTM_i0;
    
                        invSTM_j0 = [invStateSTM_j0,                   -invStateSTM_j0 * paramSTM_j0;
                                     zeros(nEstimatedParams, nStates),  eye(nEstimatedParams)];
                    else
                        invSTM_j0 = invStateSTM_j0;
                    end
                end
    
                output.iterationData{ii}.measResidualHistory = measResidualHistory;
    
                % ----------------------------------------------------------------------------------
                
                % Compute postfit state deviation and covariance at initial time
                % (i.e., map postfit state deviation and covariance at final measurement time to initial time)
                invSTM_i0 = invSTM_j0;
    
                postAugStateDelta_0 = invSTM_i0 * postAugStateDelta_i;
                postAugStateCovar_0 = invSTM_i0 * postAugStateCovar_i * invSTM_i0';
    
                if ii == 1
                    % Determine if state has converged
                    if max(abs(postAugStateDelta_0 ./ priorAugState_0)) < Settings.DEFAULT_CONVERGENCE_TOL
                        hasConverged = true;
                    end
                    
                    % Update prefit state (now becomes postfit state)
                    postAugState_0 = priorAugState_0 + postAugStateDelta_0;
                else
                    % Determine if state has converged
                    if max(abs(postAugStateDelta_0 ./ postAugState_0)) < Settings.DEFAULT_CONVERGENCE_TOL
                        hasConverged = true;
                    end
                    
                    % Update postfit state
                    postAugState_0 = postAugState_0 + postAugStateDelta_0;
                end
    
                fprintf("%i\t\t\t", ii)
                fprintf("%.4f\t", postAugState_0(:))
                fprintf("\n")
                
                output.augState_0Iterations(:, ii + 1) = postAugState_0;
                output.augStateCovar_0Iterations(:, ii + 1) = postAugStateCovar_0(:);
                
                % Shift prefit state deviation at initial time
                priorAugStateDelta_0 = priorAugStateDelta_0 - postAugStateDelta_0;
    
                % ----------------------------------------------------------------------------------
                
                % Extract postfit projectile state at initial time
                postState_0 = postAugState_0(1:nStates);
                postStateCovar_0 = postAugStateCovar_0(1:nStates, 1:nStates);
    
                output.state_0Iterations(:, ii + 1) = postState_0;
                output.stateCovar_0Iterations(:, ii + 1) = postStateCovar_0(:);
                
                % Update projectile model state (for nominal trajectory on next iteration)
                self.projectileModel.time = 0;  % See TODO
                self.projectileModel.state = postState_0;
    
                if self.includeParamSTM
                    % Extract postfit parameters
                    postParams = postAugState_0((nStates + 1):end);
                    postParamCovar = postAugStateCovar_0((nStates + 1):end, (nStates + 1):end);
    
                    output.paramIterations(:, ii + 1) = postParams;
                    output.paramCovarIterations(:, ii + 1) = postParamCovar(:);
    
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

            output.iterationData{nIterations + 1}.measResidualHistory = measResidualHistory;
            
            % Remove all unused entries in output data if converged early
            if hasConverged
                output.iterationData((nIterations + 2):end) = [];

                output.state_0Iterations(:, (nIterations + 2):end) = [];
                output.stateCovar_0Iterations(:, (nIterations + 2):end) = [];

                output.paramIterations(:, (nIterations + 2):end) = [];
                output.paramCovarIterations(:, (nIterations + 2):end) = [];

                output.augState_0Iterations(:, (nIterations + 2):end) = [];
                output.augStateCovar_0Iterations(:, (nIterations + 2):end) = [];
            end
        end
    end
end

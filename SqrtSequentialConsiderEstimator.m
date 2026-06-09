classdef SqrtSequentialConsiderEstimator < Estimator
    % TODO: Currently assumes init projectile time = estimate time epoch. Need to propagate if not

    % TODO: If P0 semi-definite (i.e., has zeros on diagonal), should really take square root by
    % doing chol() on non-zero block diagonal parts, then reassembling back with zero diagonal parts
    % Ex: P0 = [0,            ->  W0 = [0,
    %              1, 0.5                  1,
    %              0.5, 1,                 0.5, 0.8660,
    %                      0]                           0]
    % Don't need to worry about cases with zero variance + non-zero cross-covariances since states
    % with zero variance are fully known (no uncertainty), thus corresponding cross-covariances are
    % also zero by definition
    % Current band-aid with .^ 0.5 works only if P0 is diagonal (i.e., zero cross-covariances)

    methods
        % Constructor ==============================================================================

        function self = SqrtSequentialConsiderEstimator(varargin)
            self = self@Estimator(varargin{:});
        end


        % Solve method =============================================================================

        function output = solve(self, measHistory)
            % Remove any measurements from sensors not included in sensor model array
            measHistory = measHistory(:, ismember(measHistory(1, :), self.sensorModelIDs));

            nSamples = size(measHistory, 2);
            
            measTimeHistory = measHistory(2, :);   % Get measurement time history
            finalMeasTime = measTimeHistory(end);  % Get time of final measurement

            nSensorModels = length(self.sensorModelArray);

            % --------------------------------------------------------------------------------------

            nStates = self.projectileModel.nStates;

            nEstimatedProjectileParams = self.projectileModel.nEstimatedParams;
            nEstimatedPlanetParams     = self.planetModel.nEstimatedParams;
            nEstimatedParams           = nEstimatedProjectileParams + nEstimatedPlanetParams;

            nConsideredProjectileParams = self.projectileModel.nConsideredParams;
            nConsideredPlanetParams     = self.planetModel.nConsideredParams;
            nConsideredDynamicsParams   = nConsideredProjectileParams + nConsideredPlanetParams;
            nConsideredSensorParams     = 0;
            for i = 1:nSensorModels
                nConsideredSensorParams = nConsideredSensorParams + self.sensorModelArray{i}.nConsideredParams;
            end
            nConsideredParams           = nConsideredDynamicsParams + nConsideredSensorParams;

            self.includeParamSTM                         = logical(nEstimatedParams);
            self.projectileModelDynamics.includeParamSTM = logical(nEstimatedParams);

            self.includeConsideredDynamicsParamSTM                         = logical(nConsideredDynamicsParams);
            self.projectileModelDynamics.includeConsideredDynamicsParamSTM = logical(nConsideredDynamicsParams);

            nAugStates         = nStates + nEstimatedParams;

            % Build prefit vectors and covariances at initial time
            priorState_0      = self.projectileModel.state;
            priorStateCovar_0 = self.projectileModel.stateCovar;

            priorParams     = [self.projectileModel.estimatedParams; self.planetModel.estimatedParams];
            priorParamCovar = blkdiag(self.projectileModel.estimatedParamCovar, self.planetModel.estimatedParamCovar);

            priorConsideredParamCovar = blkdiag(self.projectileModel.consideredParamCovar, self.planetModel.consideredParamCovar);
            for i = 1:nSensorModels
                priorConsideredParamCovar = blkdiag(priorConsideredParamCovar, self.sensorModelArray{i}.consideredParamCovar);
            end

            priorAugState_0      = [priorState_0; priorParams];
            priorAugStateCovar_0 = blkdiag(priorStateCovar_0, priorParamCovar);
            try
                priorAugStateCovarSqrt_0 = chol(priorAugStateCovar_0)';  % P0 is positive definite

                priorConsideredParamCovarSqrt = chol(priorConsideredParamCovar);

            catch
                priorAugStateCovarSqrt_0 = priorAugStateCovar_0 .^ 0.5;  % P0 is positive semi-definite. See TODO

                priorConsideredParamCovarSqrt = priorConsideredParamCovar .^ 0.5;
            end
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
                if self.includeConsideredDynamicsParamSTM
                    [nomTimeHistory, nomStateHistory, stateSTMHistory, paramSTMHistory, consideredDynamicsParamSTMHistory] = ...
                        self.propagator.propagateWithConsideredSTM(finalMeasTime);
                else
                    [nomTimeHistory, nomStateHistory, stateSTMHistory, paramSTMHistory] = ...
                        self.propagator.propagateWithSTM(finalMeasTime);
                end
                
                output.perIterationData{ii}.nomTimeHistory  = nomTimeHistory;
                output.perIterationData{ii}.nomStateHistory = nomStateHistory;
                
                % Resample nominal trajectory at measurement times (should be exact)
                nomStateHistory = Utils.resampleStateHistory(nomTimeHistory, nomStateHistory, measTimeHistory);
                stateSTMHistory = Utils.resampleStateHistory(nomTimeHistory, stateSTMHistory, measTimeHistory);
                if self.includeParamSTM
                    paramSTMHistory = Utils.resampleStateHistory(nomTimeHistory, paramSTMHistory, measTimeHistory);
                end
                if self.includeConsideredDynamicsParamSTM
                    consideredDynamicsParamSTMHistory = Utils.resampleStateHistory(nomTimeHistory, consideredDynamicsParamSTMHistory, measTimeHistory);
                end
                
                % ----------------------------------------------------------------------------------
                
                % Note: t_i = current measurement time
                %       t_j = t_(i-1) = previous measurement time
                
                % Initialize postfit state deviation and covariance sqrt at previous measurement time (i.e., initial time here)
                postAugStateDelta_j     = priorAugStateDelta_0;
                postAugStateCovarSqrt_j = priorAugStateCovarSqrt_0;
                
                postConsideredParamSensitivity_j = zeros(nAugStates, nConsideredParams);
                
                % Initialize STMs at previous measurement time (i.e., initial time here)
                invStateSTM_j0 = eye(nStates);
                if self.includeParamSTM
                    paramSTM_j0 = zeros(nStates, nEstimatedParams);
    
                    invSTM_j0 = [invStateSTM_j0,                   -invStateSTM_j0 * paramSTM_j0;
                                 zeros(nEstimatedParams, nStates),  eye(nEstimatedParams)];
                else
                    invSTM_j0 = invStateSTM_j0;
                end

                consideredSTM_j0 = zeros(nAugStates, nConsideredParams);
    
                % Initialize measurement residual history
                measResidualHistory       = zeros(size(measHistory, 1), nSamples);
                measResidualHistory(1, :) = measHistory(1, :);
                measResidualHistory(2, :) = measHistory(2, :);

                postAugStateCovarHistory              = zeros(nAugStates ^ 2, nSamples);
                postConsideredAugStateCovarHistory    = zeros(nAugStates ^ 2, nSamples);
                postConsideredParamCrossCovarHistory  = zeros(nAugStates * nConsideredParams, nSamples);
                postConsideredParamSensitivityHistory = zeros(nAugStates * nConsideredParams, nSamples);
    
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
                    if self.includeConsideredDynamicsParamSTM
                        consideredDynamicsParamSTM_i0 = consideredDynamicsParamSTMHistory(:, i);
                        consideredDynamicsParamSTM_i0 = reshape(consideredDynamicsParamSTM_i0, [nStates, nConsideredDynamicsParams]);
                    end
    
                    if self.includeParamSTM
                        STM_i0 = [stateSTM_i0,                      paramSTM_i0;
                                  zeros(nEstimatedParams, nStates), eye(nEstimatedParams)];
                    else
                        STM_i0 = stateSTM_i0;
                    end
                    
                    if self.includeConsideredDynamicsParamSTM
                        consideredSTM_i0 = [consideredDynamicsParamSTM_i0,                      zeros(nStates, nConsideredSensorParams);
                                            zeros(nEstimatedParams, nConsideredDynamicsParams), zeros(nEstimatedParams, nConsideredSensorParams)];
                    else
                        consideredSTM_i0 = zeros(nAugStates, nConsideredParams);
                    end
                    
                    % Compute step STM from previous measurement time to current measurement time
                    STM_ij = STM_i0 * invSTM_j0;

                    consideredSTM_ij = consideredSTM_i0 - STM_ij * consideredSTM_j0;
    
                    % Propagate prefit state deviation and covariance sqrt to current measurement time
                    priorAugStateDelta_i = STM_ij * postAugStateDelta_j;
                    priorAugStateCovarSqrt_i = qr((STM_ij * postAugStateCovarSqrt_j)')';

                    priorConsideredParamSensitivity_i = STM_ij * postConsideredParamSensitivity_j + consideredSTM_ij;

                    % priorConsideredAugStateCovarSqrt_i   = qr([priorAugStateCovarSqrt_i, priorConsideredParamSensitivity_i * priorConsideredParamCovarSqrt]')';
                    % priorConsideredAugStateCovarSqrt_i   = priorConsideredAugStateCovarSqrt_i(1:nAugStates, 1:nAugStates);
                    % priorConsideredParamCrossCovarSqrt_i = priorConsideredParamSensitivity_i * priorConsideredParamCovarSqrt;
                    
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

                    consideredParamH_i = sensorModel_i.computeConsideredParamJacobian(nomState_i);
                    
                    % Get measurement noise covariance sqrt
                    measNoiseCovarSqrt_i = sensorModel_i.measNoiseCovarSqrt;

                    % Construct joint innovation + prefit state covariance sqrt
                    priorJointCovarSqrt_i = [measNoiseCovarSqrt_i,     H_i * priorAugStateCovarSqrt_i;
                                             zeros(nAugStates, nMeas), priorAugStateCovarSqrt_i];
                    
                    % Compute joint innovation + postfit state covariance sqrt
                    postJointCovarSqrt_i = qr(priorJointCovarSqrt_i')';

                    % Extract innovation covariance sqrt
                    innovCovarSqrt_i = postJointCovarSqrt_i(1:nMeas, 1:nMeas);
                    mappedInnovCovarSqrt_i = postJointCovarSqrt_i((nMeas + 1):end, 1:nMeas);
                    
                    % Compute filter gain (i.e., Kalman gain) matrix
                    measResidualGain_i = mappedInnovCovarSqrt_i / innovCovarSqrt_i;
                    
                    % Update to postfit state deviation and covariance sqrt using current measurement residual
                    postAugStateDelta_i     = priorAugStateDelta_i + measResidualGain_i * (measResidual_i - H_i * priorAugStateDelta_i);
                    postAugStateCovarSqrt_i = postJointCovarSqrt_i((nMeas + 1):end, (nMeas + 1):end);
                    
                    postConsideredParamSensitivity_i = (eye(nAugStates) - measResidualGain_i * H_i) * priorConsideredParamSensitivity_i - measResidualGain_i * consideredParamH_i;

                    postConsideredAugStateCovarSqrt_i   = qr([postAugStateCovarSqrt_i, postConsideredParamSensitivity_i * priorConsideredParamCovarSqrt]')';
                    postConsideredAugStateCovarSqrt_i   = postConsideredAugStateCovarSqrt_i(1:nAugStates, 1:nAugStates);
                    postConsideredParamCrossCovarSqrt_i = postConsideredParamSensitivity_i * priorConsideredParamCovarSqrt;

                    postAugStateCovarHistory(:, i)              = reshape(postAugStateCovarSqrt_i * postAugStateCovarSqrt_i', [nAugStates ^ 2, 1]);
                    postConsideredAugStateCovarHistory(:, i)    = reshape(postConsideredAugStateCovarSqrt_i * postConsideredAugStateCovarSqrt_i', [nAugStates ^ 2, 1]);
                    postConsideredParamCrossCovarHistory(:, i)  = reshape(postConsideredParamCrossCovarSqrt_i * priorConsideredParamCovarSqrt', [nAugStates * nConsideredParams, 1]);
                    postConsideredParamSensitivityHistory(:, i) = reshape(postConsideredParamSensitivity_i, [nAugStates * nConsideredParams, 1]);

                    % Store results (current measurement time now becomes previous measurement time)
                    postAugStateDelta_j = postAugStateDelta_i;
                    postAugStateCovarSqrt_j = postAugStateCovarSqrt_i;

                    postConsideredParamSensitivity_j = postConsideredParamSensitivity_i;
                    
                    invStateSTM_j0 = inv(stateSTM_i0);
                    if self.includeParamSTM
                        paramSTM_j0 = paramSTM_i0;
    
                        invSTM_j0 = [invStateSTM_j0,                   -invStateSTM_j0 * paramSTM_j0;
                                     zeros(nEstimatedParams, nStates),  eye(nEstimatedParams)];
                    else
                        invSTM_j0 = invStateSTM_j0;
                    end

                    consideredSTM_j0 = consideredSTM_i0;
                end
    
                output.perIterationData{ii}.measTimeHistory = measTimeHistory;
                output.perIterationData{ii}.measResidualHistory = measResidualHistory;

                output.perIterationData{ii + 1}.postAugStateCovarHistory              = postAugStateCovarHistory;
                output.perIterationData{ii + 1}.postConsideredAugStateCovarHistory    = postConsideredAugStateCovarHistory;
                output.perIterationData{ii + 1}.postConsideredParamCrossCovarHistory  = postConsideredParamCrossCovarHistory;
                output.perIterationData{ii + 1}.postConsideredParamSensitivityHistory = postConsideredParamSensitivityHistory;
    
                % ----------------------------------------------------------------------------------
                
                % Compute postfit state deviation and covariance at initial time
                % (i.e., map postfit state deviation and covariance sqrt at final measurement time to initial time)
                invSTM_i0 = invSTM_j0;
    
                postAugStateDelta_0 = invSTM_i0 * postAugStateDelta_i;
                postAugStateCovarSqrt_0 = qr((invSTM_i0 * postAugStateCovarSqrt_i)')';
                postAugStateCovar_0 = postAugStateCovarSqrt_0 * postAugStateCovarSqrt_0';
    
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

            output.perIterationData{nIterations + 1}.measTimeHistory = measTimeHistory;
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

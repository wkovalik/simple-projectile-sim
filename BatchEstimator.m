classdef BatchEstimator < handle
    properties
        projectileModelDynamics

        projectileModel
        earthModel

        integrator
        propagator

        estimatedParams
        estimatedProjectileParams
        estimatedEarthParams

        sensorModelList
        sensorModelMap

        stateInterpolantList
        measurementResidualRecordList
    end

    properties (Dependent)
        N_STATES
        N_ESTIMATED_PARAMS
        N_ESTIMATED_PROJECTILE_PARAMS
        N_ESTIMATED_EARTH_PARAMS
        N_SENSOR_MODELS
    end

    methods
        % Constructor ==============================================================================
        function self = BatchEstimator(projectileModelDynamics, integrator, sensorModelList)
            % Set handle for estimator projectile model dynamics
            self.projectileModelDynamics = projectileModelDynamics;
            
            % Set handle for estimator projectile and planet models (from dynamics model)
            self.projectileModel = self.projectileModelDynamics.projectile;
            self.earthModel = self.projectileModelDynamics.earth;
            
            % Set handle for estimator integrator
            self.integrator = integrator;
            
            % Create estimator propagator
            self.propagator = Propagator(self.projectileModelDynamics, self.integrator);
            self.propagator.selectPropagator("stateAndSTM");  % Option to propagate state, state STM, and parameter STM (required for estimation)
            
            % Store list of sensor models to be used during estimation
            self.sensorModelList = sensorModelList;
            
            % Process list of sensor models
            self.sensorModelMap = dictionary();
            for i = 1:self.N_SENSOR_MODELS
                % Pass handle for projectile model to each sensor model
                self.sensorModelList{i}.projectileModel = self.projectileModel;
                
                % Add sensor model to map of sensor model IDs
                self.sensorModelMap = insert( ...
                    self.sensorModelMap, self.sensorModelList{i}.sensorID, self.sensorModelList(i) ...
                );
            end

            % Preallocate lists to store data on each estimator iteration
            self.stateInterpolantList = cell(1, Constants.MAX_ESTIMATOR_ITERATIONS + 1);           % List to store estimated state trajectory interpolants
            self.measurementResidualRecordList = cell(1, Constants.MAX_ESTIMATOR_ITERATIONS + 1);  % List to store measurement residuals
        end
        

        % Solve method =============================================================================
        function solve(self, estimateTime, sensorSampleRecord)
            % Setup and preprocessing --------------------------------------------------------------

            % Remove measurements from any sensors which are not included in sensor model list
            sensorSampleRecord = self.removeUnmodeledSensorsFromRecord(sensorSampleRecord);

            % Create aggregate measurement history with measurements sorted in chronological order
            [chronoSensorSampleRecord, nTotalSamples, finalSampleTime] = ...
                self.createChronoSensorSampleRecord(sensorSampleRecord);

            % Find all projectile and planet parameters to be estimated
            self.findEstimatedParams();
            
            % Preallocate list to store measurement residuals on each iteration
            self.initMeasurementResidualRecord(sensorSampleRecord);


            % Estimator state initialization -------------------------------------------------------

            % Propagate projectile state and covariance to desired estimate time
            self.projectileModel.time = estimateTime;  % TODO: Currently assumes estimateTime = initSimTime. Else, would need to propagate projectile state and covariance to estimateTime
            
            % Store prefit projectile state and covariance
            priorState = self.projectileModel.state.vector;
            priorStateCovar = self.projectileModel.state.covar;
            
            % Store prefit projectile and planet model parameters and covariance
            [priorParamState, priorParamCovar] = self.getEstimatedParamState();

            % Create augmented prefit state and covariance
            priorAugState = [priorState; priorParamState];
            priorAugStateCovar = blkdiag(priorStateCovar, priorParamCovar);

            % Initialize prefit state deviation
            priorAugStateDelta = zeros(self.N_STATES + self.N_ESTIMATED_PARAMS, 1);
            

            % Main estimation loop -----------------------------------------------------------------

            for ii = 1:Constants.MAX_ESTIMATOR_ITERATIONS  % TODO: Implement covergence break out
                % Batch time update  ---------------------------------------------------------------

                % Propagate nominal projectile state and STMs
                self.propagator.propagate(finalSampleTime);
                self.stateInterpolantList{ii} = self.propagator.stateInterpolant;  % Store nominal projectile trajectory
    
                % Batch measurement update  --------------------------------------------------------

                % Initialize postfit information matrix and normal vector
                postAugStateInfor = inv(priorAugStateCovar);
                postNormalVector = priorAugStateCovar \ priorAugStateDelta;
    
                for i = 1:nTotalSamples
                    % Get next sensor measurement time and measurement
                    sampleTime = chronoSensorSampleRecord{i}.sampleTime;
                    sampledMeasurement = chronoSensorSampleRecord{i}.sampledMeasurement;
    
                    % Get nominal projectile state and STMs at time of measurement
                    sampleTimeState = self.propagator.getState(sampleTime);
                    [sampleTimeStateSTM, sampleTimeParamSTM] = self.propagator.getSTM(sampleTime);
    
                    % Get model for next sensor measurement
                    sensorID = chronoSensorSampleRecord{i}.sensorID;
                    sensorModel = self.sensorModelMap{sensorID};
                    
                    % Compute measurement according to sensor model at nominal state
                    computedMeasurement = sensorModel.computeMeasurement(sampleTimeState);
                    
                    % Compute measurement residual
                    measurementResidual = sampledMeasurement - computedMeasurement;
                    self.addResidualToRecord(sensorID, sampleTime, measurementResidual, ii);  % Store measurement residual
                    
                    % Compute measurement Jacobians
                    stateH = sensorModel.computeStateJacobian(sampleTimeState);
                    mappedStateH = stateH * sampleTimeStateSTM;

                    if ~isempty(self.estimatedParams)
                        paramH = sensorModel.computeParamJacobian(sampleTimeState);
                        mappedParamH = stateH * sampleTimeParamSTM + paramH;
                    else
                        mappedParamH = [];
                    end
    
                    mappedH = [mappedStateH, mappedParamH];
    
                    % Compute next contribution to postfit information matrix and normal vector from current measurement
                    nextAugStateInfor = mappedH' * sensorModel.invNoiseCovar * mappedH;
                    nextNormalVector = mappedH' * sensorModel.invNoiseCovar * measurementResidual;
                    
                    % Add contribution to postfit information matrix and normal vector
                    postAugStateInfor = postAugStateInfor + nextAugStateInfor;
                    postNormalVector = postNormalVector + nextNormalVector;
                end
                
                % Compute postfit state deviation and covariance
                postAugStateDelta = postAugStateInfor \ postNormalVector;
                postAugStateCovar = inv(postAugStateInfor);
                
                % Compute postfit state (and update on subsequent iterations)
                if ii == 1
                    postAugState = priorAugState + postAugStateDelta;
                else
                    postAugState = postAugState + postAugStateDelta;
                end
                
                % Shift prefit state deviation to current postfit state
                priorAugStateDelta = priorAugStateDelta - postAugStateDelta;

                % Update projectile using postfit state and covariance (becomes new nominal state)
                postState = postAugState(1:self.N_STATES);
                postStateCovar = postAugStateCovar(1:self.N_STATES, 1:self.N_STATES);
                
                self.projectileModel.updateState(estimateTime, postState, postStateCovar);
                
                % Update projectile and planet models using postfit parameters and covariance
                if ~isempty(self.estimatedParams)
                    postParamState = postAugState(self.N_STATES + (1:self.N_ESTIMATED_PARAMS));
                    postParamCovar = postAugStateCovar(self.N_STATES + (1:self.N_ESTIMATED_PARAMS), self.N_STATES + (1:self.N_ESTIMATED_PARAMS));

                    self.updateEstimatedParams(postParamState, postParamCovar);
                end
            end

            % --------------------------------------------------------------------------------------

            % Postfit processing -------------------------------------------------------------------

            % Propagate postfit projectile state
            self.propagator.propagate(finalSampleTime);
            self.stateInterpolantList{ii + 1} = self.propagator.stateInterpolant;  % Store postfit projectile trajectory
            
            % Compute postfit measurement residuals
            for i = 1:nTotalSamples
                % Get next sensor measurement time and measurement
                sampleTime = chronoSensorSampleRecord{i}.sampleTime;
                sampledMeasurement = chronoSensorSampleRecord{i}.sampledMeasurement;
                
                % Get postfit projectile state at time of measurement
                sampleTimeState = self.propagator.getState(sampleTime);
                
                % Get model for next sensor measurement 
                sensorID = chronoSensorSampleRecord{i}.sensorID;
                sensorModel = self.sensorModelMap{sensorID};
                
                % Compute measurement according to sensor model at postfit state
                computedMeasurement = sensorModel.computeMeasurement(sampleTimeState);
                
                % Compute measurement residual
                measurementResidual = sampledMeasurement - computedMeasurement;
                self.addResidualToRecord(sensorID, sampleTime, measurementResidual, ii + 1);  % Store measurement residual
            end
        end


        % Helper methods ===========================================================================
        function findEstimatedParams(self)
            % Find all projectile parameters set to be estimated
            self.estimatedProjectileParams = {};
            
            projectileParamNames = fieldnames(self.projectileModel.params);
            for i = 1:length(projectileParamNames)
                paramName = projectileParamNames{i};
            
                if self.projectileModel.params.(paramName).isEstimated
                    self.estimatedProjectileParams(end + 1) = { paramName };
                end
            end

            self.estimatedProjectileParams = string(self.estimatedProjectileParams);
            
            % Find all planet parameters set to be estimated
            self.estimatedEarthParams = {};

            earthParamNames = fieldnames(self.earthModel.params);
            for i = 1:length(earthParamNames)
                paramName = earthParamNames{i};
            
                if self.earthModel.params.(paramName).isEstimated
                    self.estimatedEarthParams(end + 1) = { paramName };
                end
            end
            
            self.estimatedEarthParams = string(self.estimatedEarthParams);

            % Concatenate and store list of parameters to be estimated
            self.estimatedParams = [self.estimatedProjectileParams, self.estimatedEarthParams];            
            self.propagator.estimatedParams = self.estimatedParams;

            % Pass list to projectile model dynamics
            self.projectileModelDynamics.estimatedParams = self.estimatedParams;
            
            % Pass list to sensor models
            for i = 1:length(self.sensorModelList)
                self.sensorModelList{i}.estimatedParams = self.estimatedParams;
            end
        end

        function [paramState, paramCovar] = getEstimatedParamState(self)
            % Build vector of estimated parameter values, as well as covariance matrix

            if ~isempty(self.estimatedParams)
                paramState = zeros(self.N_ESTIMATED_PARAMS, 1);
                paramCovar = zeros(self.N_ESTIMATED_PARAMS);
    
                % Get estimated projectile parameters and covariances
                for i = 1:self.N_ESTIMATED_PROJECTILE_PARAMS
                    paramName = self.estimatedProjectileParams(i);
    
                    paramState(i) = self.projectileModel.params.(paramName).value;
                    paramCovar(i, i) = self.projectileModel.params.(paramName).covar;
                end
                
                % Get estimated planet parameters and covariances
                for i = 1:self.N_ESTIMATED_EARTH_PARAMS
                    paramName = self.estimatedEarthParams(i);
    
                    paramState(i + self.N_ESTIMATED_PROJECTILE_PARAMS) = self.earthModel.params.(paramName).value;
                    paramCovar(i + self.N_ESTIMATED_PROJECTILE_PARAMS, i + self.N_ESTIMATED_PROJECTILE_PARAMS) = self.earthModel.params.(paramName).covar;
                end
                
            else
                paramState = [];
                paramCovar = [];
            end
        end

        function updateEstimatedParams(self, paramState, paramCovar)
            % Use postfit estimated parameter vector to update corresponding values in projectile
            % and planet models

            % Update projectile model parameters and covariances
            for i = 1:self.N_ESTIMATED_PROJECTILE_PARAMS
                paramName = self.estimatedProjectileParams(i);

                self.projectileModel.params.(paramName).value = paramState(i);
                self.projectileModel.params.(paramName).covar = paramCovar(i, i);
            end

            % Update planet model parameters and covariances
            for i = 1:self.N_ESTIMATED_EARTH_PARAMS
                paramName = self.estimatedEarthParams(i);

                self.earthModel.params.(paramName).value = paramState(i + self.N_ESTIMATED_PROJECTILE_PARAMS);
                self.earthModel.params.(paramName).covar = paramCovar(i + self.N_ESTIMATED_PROJECTILE_PARAMS, i + self.N_ESTIMATED_PROJECTILE_PARAMS);
            end
        end


        % Sensor methods ===========================================================================
        function addSensorModel(self, sensorModel)
            % Append new sensor model to list
            self.sensorModelList(end + 1) = { sensorModel };
                
            % Pass projectile model to new sensor model
            self.sensorModelList{end}.projectileModel = self.projectileModel;
            
            % Add new sensor model to map of sensor model IDs
            self.sensorModelMap = insert( ...
                self.sensorModelMap, self.sensorModelList{end}.sensorID, self.sensorModelList(end) ...
            );
        end

        function sensorSampleRecord = removeUnmodeledSensorsFromRecord(self, sensorSampleRecord)
            % Remove measurements from any sensors which are not included in sensor model list

            sensorIDs = keys(sensorSampleRecord);

            for i = 1:length(sensorIDs)
                sensorID = sensorIDs(i);
                
                isUnmodeledSensor = true;
                for j = 1:self.N_SENSOR_MODELS
                    sensorModelID = self.sensorModelList{j}.sensorID;

                    if sensorID == sensorModelID
                        isUnmodeledSensor = false;
                    end
                end

                if isUnmodeledSensor
                    sensorSampleRecord = remove(sensorSampleRecord, sensorID);
                end
            end
        end

        function [chronoSensorSampleRecord, nTotalSamples, finalSampleTime] = createChronoSensorSampleRecord( ...
                ~, sensorSampleRecord ...
        )
            % Create aggregate measurement history with measurements sorted in chronological order

            sensorIDs = keys(sensorSampleRecord);
        
            % Find total number of measurements from all sensors
            nTotalSamples = 0;
            for i = 1:length(sensorIDs)
                sensorID = sensorIDs(i);
                nSamples = sensorSampleRecord(sensorID).nSamples;
        
                nTotalSamples = nTotalSamples + nSamples;
            end
        
            % Build concatenated, unsorted list of measurement histories from all sensors
            concatSensorSampleRecord = cell(1, nTotalSamples);
            sampleCount = 0;
            for i = 1:length(sensorIDs)
                sensorID = sensorIDs(i);
                nSamples = sensorSampleRecord(sensorID).nSamples;
        
                for j = 1:nSamples
                    sampleCount = sampleCount + 1;
                    
                    concatSensorSampleRecord{sampleCount}.sensorID = sensorID;
                    concatSensorSampleRecord{sampleCount}.sampleTime = sensorSampleRecord(sensorID).timeHistory(j);
                    concatSensorSampleRecord{sampleCount}.sampledMeasurement = sensorSampleRecord(sensorID).sampleHistory(:, j);
                end
            end
            
            % Determine order to sort concatenated list in chronological order
            sampleTimes = zeros(1, nTotalSamples);
            sensorIDs = zeros(1, nTotalSamples);
            for i = 1:nTotalSamples
                sampleTimes(i) = concatSensorSampleRecord{i}.sampleTime;
                sensorIDs(i) = concatSensorSampleRecord{i}.sensorID;
            end
        
            [~, sortIndices] = sortrows([sampleTimes; sensorIDs]', [1, 2]);  % Sort by measurement time, then by sensor ID
        
            % Sort concatenated list in chronological order
            chronoSensorSampleRecord = cell(1, nTotalSamples);
            for i = 1:nTotalSamples
                chronoSensorSampleRecord{i} = concatSensorSampleRecord{sortIndices(i)};
            end

            finalSampleTime = chronoSensorSampleRecord{end}.sampleTime;
            % finalSampleTime = self.getFinalSampleTime(sensorSampleRecord);
        end

        function finalSampleTime = getFinalSampleTime(~, sensorSampleRecord)
            % Get time of final sensor measurement

            sensorIDs = keys(sensorSampleRecord);

            finalSampleTime = 0;
            for i = 1:length(sensorIDs)
                sensorID = sensorIDs(i);
                sensorFinalSampleTime = sensorSampleRecord(sensorID).timeHistory(end);
            
                if sensorFinalSampleTime > finalSampleTime
                    finalSampleTime = sensorFinalSampleTime;
                end
            end
        end

        function initMeasurementResidualRecord(self, sensorSampleRecord)
            % Preallocate list to store measurement residuals on each iteration

            for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
                self.measurementResidualRecordList{ii} = dictionary();
    
                for i = 1:self.N_SENSOR_MODELS
                    sensorModel = self.sensorModelList{i};
                    
                    nSamples = sensorSampleRecord(sensorModel.sensorID).nSamples;
                    
                    initRecord.nSamples = 0;
                    initRecord.timeHistory = zeros(1, nSamples);
                    initRecord.residualHistory = zeros(sensorModel.N_MEASUREMENTS, nSamples);
    
                    self.measurementResidualRecordList{ii} = insert( ...
                        self.measurementResidualRecordList{ii}, sensorModel.sensorID, initRecord ...
                    );
                end
            end
        end

        function addResidualToRecord(self, sensorID, sampleTime, measurementResidual, nIteration)
            % Add measurement residual to history record 

            nSamples = self.measurementResidualRecordList{nIteration}(sensorID).nSamples;
            nSamples = nSamples + 1;
            
            self.measurementResidualRecordList{nIteration}(sensorID).nSamples = nSamples;
            self.measurementResidualRecordList{nIteration}(sensorID).timeHistory(nSamples) = sampleTime;
            self.measurementResidualRecordList{nIteration}(sensorID).residualHistory(:, nSamples) = measurementResidual;
        end


        % Getters ==================================================================================
        function N_STATES = get.N_STATES(self)
            N_STATES = self.projectileModel.state.N_STATES;
        end

        function N_ESTIMATED_PARAMS = get.N_ESTIMATED_PARAMS(self)
            N_ESTIMATED_PARAMS = length(self.estimatedParams);
        end

        function N_ESTIMATED_PROJECTILE_PARAMS = get.N_ESTIMATED_PROJECTILE_PARAMS(self)
            N_ESTIMATED_PROJECTILE_PARAMS = length(self.estimatedProjectileParams);
        end

        function N_ESTIMATED_EARTH_PARAMS = get.N_ESTIMATED_EARTH_PARAMS(self)
            N_ESTIMATED_EARTH_PARAMS = length(self.estimatedEarthParams);
        end

        function N_SENSOR_MODELS = get.N_SENSOR_MODELS(self)
            N_SENSOR_MODELS = length(self.sensorModelList);
        end
    end
end

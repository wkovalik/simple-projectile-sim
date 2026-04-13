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
            self.projectileModelDynamics = projectileModelDynamics;

            self.projectileModel = self.projectileModelDynamics.projectile;
            self.earthModel = self.projectileModelDynamics.earth;

            self.integrator = integrator;

            self.propagator = Propagator(self.projectileModelDynamics, self.integrator);
            self.propagator.selectPropagator("stateAndSTM");

            self.sensorModelList = sensorModelList;
            self.sensorModelMap = dictionary();

            for i = 1:self.N_SENSOR_MODELS
                self.sensorModelList{i}.projectileModel = self.projectileModel;
                
                self.sensorModelMap = insert( ...
                    self.sensorModelMap, self.sensorModelList{i}.sensorID, self.sensorModelList(i) ...
                );
            end

            self.stateInterpolantList = cell(1, Constants.MAX_ESTIMATOR_ITERATIONS + 1);
            self.measurementResidualRecordList = cell(1, Constants.MAX_ESTIMATOR_ITERATIONS + 1);
        end
        
        % Solve method =============================================================================
        function solve(self, estimateTime, sensorSampleRecord)
            sensorSampleRecord = self.removeUnmodeledSensorsFromRecord(sensorSampleRecord);
            [chronoSensorSampleRecord, nTotalSamples, finalSampleTime] = ...
                self.createChronoSensorSampleRecord(sensorSampleRecord);

            self.findEstimatedParams();

            self.projectileModel.time = estimateTime;  % TODO: Propagate to estTime if not starting at estTime
            
            priorState = self.projectileModel.state.vector;
            priorStateCovar = self.projectileModel.state.covar;

            [priorParamState, priorParamCovar] = self.getEstimatedParamState();

            priorAugState = [priorState; priorParamState];

            priorAugStateCovar = blkdiag(priorStateCovar, priorParamCovar);
            priorAugStateDelta = zeros(self.N_STATES + self.N_ESTIMATED_PARAMS, 1);
            
            for ii = 1:Constants.MAX_ESTIMATOR_ITERATIONS  % TODO: Implement covergence break out
                self.initMeasurementResidualRecord(sensorSampleRecord, ii);

                % Propagate nominal trajectory and STM
                self.propagator.propagate(finalSampleTime);

                self.stateInterpolantList{ii} = self.propagator.stateInterpolant;
    
                % Process measurements
                postAugStateInfor = inv(priorAugStateCovar);
                postNormalVector = priorAugStateCovar \ priorAugStateDelta;
    
                for i = 1:nTotalSamples
                    sensorID = chronoSensorSampleRecord{i}.sensorID;
                    sampleTime = chronoSensorSampleRecord{i}.sampleTime;
                    sampledMeasurement = chronoSensorSampleRecord{i}.sampledMeasurement;
    
                    % Get nominal state and STM
                    sampleTimeState = self.propagator.getState(sampleTime);
                    [sampleTimeStateSTM, sampleTimeParamSTM] = self.propagator.getSTM(sampleTime);
    
                    % Process next measurement
                    sensorModel = self.sensorModelMap{sensorID};
    
                    computedMeasurement = sensorModel.computeMeasurement(sampleTimeState);

                    measurementResidual = sampledMeasurement - computedMeasurement;
                    self.addResidualToRecord(sensorID, sampleTime, measurementResidual, ii);
                    
                    stateH = sensorModel.computeStateJacobian(sampleTimeState);
                    mappedStateH = stateH * sampleTimeStateSTM;

                    if ~isempty(self.estimatedParams)
                        paramH = sensorModel.computeParamJacobian(sampleTimeState);
                        mappedParamH = stateH * sampleTimeParamSTM + paramH;
                    else
                        mappedParamH = [];
                    end
    
                    mappedH = [mappedStateH, mappedParamH];
    
                    addAugStateInfor = mappedH' * sensorModel.invNoiseCovar * mappedH;
                    addNormalVector = mappedH' * sensorModel.invNoiseCovar * measurementResidual;
    
                    postAugStateInfor = postAugStateInfor + addAugStateInfor;
                    postNormalVector = postNormalVector + addNormalVector;
                end
    
                postAugStateDelta = postAugStateInfor \ postNormalVector;
                postAugStateCovar = inv(postAugStateInfor);
                
                if ii == 1
                    postAugState = priorAugState + postAugStateDelta;
                else
                    postAugState = postAugState + postAugStateDelta;
                end
                priorAugStateDelta = priorAugStateDelta - postAugStateDelta;


                postState = postAugState(1:self.N_STATES);
                postStateCovar = postAugStateCovar(1:self.N_STATES, 1:self.N_STATES);
                
                self.projectileModel.updateState(estimateTime, postState, postStateCovar);
                

                if ~isempty(self.estimatedParams)
                    postParamState = postAugState(self.N_STATES + (1:self.N_ESTIMATED_PARAMS));
                    postParamCovar = postAugStateCovar(self.N_STATES + (1:self.N_ESTIMATED_PARAMS), self.N_STATES + (1:self.N_ESTIMATED_PARAMS));

                    self.updateEstimatedParams(postParamState, postParamCovar);
                end
            end

            % Postfit trajectory -------------------------------------------------------------------
            self.propagator.propagate(finalSampleTime);

            self.stateInterpolantList{ii + 1} = self.propagator.stateInterpolant;

            self.initMeasurementResidualRecord(sensorSampleRecord, ii + 1);

            for i = 1:nTotalSamples
                sensorID = chronoSensorSampleRecord{i}.sensorID;
                sampleTime = chronoSensorSampleRecord{i}.sampleTime;
                sampledMeasurement = chronoSensorSampleRecord{i}.sampledMeasurement;

                sampleTimeState = self.propagator.getState(sampleTime);

                sensorModel = self.sensorModelMap{sensorID};

                computedMeasurement = sensorModel.computeMeasurement(sampleTimeState);

                measurementResidual = sampledMeasurement - computedMeasurement;
                self.addResidualToRecord(sensorID, sampleTime, measurementResidual, ii + 1);
            end
        end

        % Helper methods ===========================================================================
        function findEstimatedParams(self)
            self.estimatedProjectileParams = {};

            projectileParamNames = fieldnames(self.projectileModel.params);
            for i = 1:length(projectileParamNames)
                paramName = projectileParamNames{i};
            
                if self.projectileModel.params.(paramName).isEstimated
                    self.estimatedProjectileParams(end + 1) = { paramName };
                end
            end

            self.estimatedProjectileParams = string(self.estimatedProjectileParams);
            
            self.estimatedEarthParams = {};

            earthParamNames = fieldnames(self.earthModel.params);
            for i = 1:length(earthParamNames)
                paramName = earthParamNames{i};
            
                if self.earthModel.params.(paramName).isEstimated
                    self.estimatedEarthParams(end + 1) = { paramName };
                end
            end
            
            self.estimatedEarthParams = string(self.estimatedEarthParams);

            self.estimatedParams = [self.estimatedProjectileParams, self.estimatedEarthParams];

            self.propagator.estimatedParams = self.estimatedParams;
            self.projectileModelDynamics.estimatedParams = self.estimatedParams;
            
            for i = 1:length(self.sensorModelList)
                self.sensorModelList{i}.estimatedParams = self.estimatedParams;
            end
        end

        function [paramState, paramCovar] = getEstimatedParamState(self)
            if ~isempty(self.estimatedParams)
                paramState = zeros(self.N_ESTIMATED_PARAMS, 1);
                paramCovar = zeros(self.N_ESTIMATED_PARAMS);
    
                for i = 1:self.N_ESTIMATED_PROJECTILE_PARAMS
                    paramName = self.estimatedProjectileParams(i);
    
                    paramState(i) = self.projectileModel.params.(paramName).value;
                    paramCovar(i, i) = self.projectileModel.params.(paramName).covar;
                end
    
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
            for i = 1:self.N_ESTIMATED_PROJECTILE_PARAMS
                paramName = self.estimatedProjectileParams(i);

                self.projectileModel.params.(paramName).value = paramState(i);
                self.projectileModel.params.(paramName).covar = paramCovar(i, i);
            end

            for i = 1:self.N_ESTIMATED_EARTH_PARAMS
                paramName = self.estimatedEarthParams(i);

                self.earthModel.params.(paramName).value = paramState(i + self.N_ESTIMATED_PROJECTILE_PARAMS);
                self.earthModel.params.(paramName).covar = paramCovar(i + self.N_ESTIMATED_PROJECTILE_PARAMS, i + self.N_ESTIMATED_PROJECTILE_PARAMS);
            end
        end

        % Sensor methods ===========================================================================
        function addSensorModel(self, sensorModel)
            self.sensorModelList(end + 1) = { sensorModel };

            self.sensorModelList{end}.projectileModel = self.projectileModel;
                
            self.sensorModelMap = insert( ...
                self.sensorModelMap, self.sensorModelList{end}.sensorID, self.sensorModelList(end) ...
            );
        end

        function sensorSampleRecord = removeUnmodeledSensorsFromRecord(self, sensorSampleRecord)
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
            sensorIDs = keys(sensorSampleRecord);
        
            nTotalSamples = 0;
            for i = 1:length(sensorIDs)
                sensorID = sensorIDs(i);
                nSamples = sensorSampleRecord(sensorID).nSamples;
        
                nTotalSamples = nTotalSamples + nSamples;
            end
        
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
        
            sampleTimes = zeros(1, nTotalSamples);
            sensorIDs = zeros(1, nTotalSamples);
            for i = 1:nTotalSamples
                sampleTimes(i) = concatSensorSampleRecord{i}.sampleTime;
                sensorIDs(i) = concatSensorSampleRecord{i}.sensorID;
            end
        
            [~, sortIndices] = sortrows([sampleTimes; sensorIDs]', [1, 2]);
        
            chronoSensorSampleRecord = cell(1, nTotalSamples);
            for i = 1:nTotalSamples
                chronoSensorSampleRecord{i} = concatSensorSampleRecord{sortIndices(i)};
            end

            finalSampleTime = chronoSensorSampleRecord{end}.sampleTime;
            % finalSampleTime = self.getFinalSampleTime(sensorSampleRecord);
        end

        function finalSampleTime = getFinalSampleTime(~, sensorSampleRecord)
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

        function initMeasurementResidualRecord(self, sensorSampleRecord, nIteration)
            self.measurementResidualRecordList{nIteration} = dictionary();

            for i = 1:self.N_SENSOR_MODELS
                sensorModel = self.sensorModelList{i};
                
                nSamples = sensorSampleRecord(sensorModel.sensorID).nSamples;
                
                initRecord.nSamples = 0;
                initRecord.timeHistory = zeros(1, nSamples);
                initRecord.residualHistory = zeros(sensorModel.N_MEASUREMENTS, nSamples);

                self.measurementResidualRecordList{nIteration} = insert( ...
                    self.measurementResidualRecordList{nIteration}, sensorModel.sensorID, initRecord ...
                );
            end
        end

        function addResidualToRecord(self, sensorID, sampleTime, measurementResidual, nIteration)
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

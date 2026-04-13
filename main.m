clear; clc; close all;

rng(0);

initSimTime = 0;
finalSimTime = 7;

estimateTime = 0;

% ==================================================================================================
% === TRUTH TRAJECTORY =============================================================================
% ==================================================================================================

% Create truth projectile --------------------------------------------------------------------------
truthProjectile = Projectile();

truthProjectile.params.m.value  = 0.145;
truthProjectile.params.S.value  = (pi / 4) * 0.075 ^ 2;
truthProjectile.params.CD.value = 0.15;

truthProjectile.time = initSimTime;
truthProjectile.state.vector = [0; 0; 0; 39.5; 0; -39.75];

% Create truth Earth -------------------------------------------------------------------------------
truthEarth = Earth();

truthEarth.params.vWindx.value = -5;

% Create truth projectile dynamics -----------------------------------------------------------------
truthProjectileDynamics = ProjectileDynamics(truthProjectile, truthEarth);

% Create range sensor ------------------------------------------------------------------------------
rangeSensor = RangeSensor();
rangeSensor.sensorID = 1;

rangeSensor.params.x.value = -1;
rangeSensor.params.y.value = 0;
rangeSensor.params.z.value = 0;

rangeSensor.samplePeriod = 0.1;

rangeSensor.noiseCovar = deg2rad(0.1) ^ 2;

% Create direction sensor --------------------------------------------------------------------------
directionSensor = DirectionSensor();
directionSensor.sensorID = 2;

directionSensor.params.x.value = -1;
directionSensor.params.y.value = 0;
directionSensor.params.z.value = 0;

directionSensor.samplePeriod = 0.01;

directionSensor.noiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);

% Create integrator --------------------------------------------------------------------------------
integrator = RK4Integrator();

integrator.stepPeriod = 0.01;

% Create truth propagator --------------------------------------------------------------------------
truthPropagator = Propagator(truthProjectileDynamics, integrator);

truthPropagator.addSensor(rangeSensor);
truthPropagator.addSensor(directionSensor);

truthPropagator.selectPropagator("stateWithSensors");

% Propagate truth trajectory -----------------------------------------------------------------------
truthPropagator.propagate(finalSimTime);

sensorSampleRecord = truthPropagator.sensorSampleRecord;



% ==================================================================================================
% === ESTIMATOR ====================================================================================
% ==================================================================================================

% Create projectile model --------------------------------------------------------------------------
projectileModel = Projectile();

projectileModel.params.m.value  = 0.145;
projectileModel.params.S.value  = (pi / 4) * 0.075 ^ 2;
projectileModel.params.CD.value = 0.20;
projectileModel.params.CD.covar = 0.1 ^ 2;

projectileModel.time = initSimTime;
projectileModel.state.vector = [0; 0; 0; 40; 0; -40];
projectileModel.state.covar = diag([0.1; 0.1; 0.1; 0.5; 0.5; 0.5] .^ 2);

% Create Earth model -------------------------------------------------------------------------------
earthModel = Earth();

earthModel.params.vWindx.value = 0;
earthModel.params.vWindx.covar = 25 ^ 2;

% Create projectile model dynamics -----------------------------------------------------------------
projectileModelDynamics = ProjectileDynamics(projectileModel, earthModel);

% Create range sensor model ------------------------------------------------------------------------
rangeSensorModel = RangeSensorModel();
rangeSensorModel.sensorID = 1;

rangeSensorModel.params.x.value = -1;
rangeSensorModel.params.y.value = 0;
rangeSensorModel.params.z.value = 0;

rangeSensorModel.noiseCovar = 0.01 ^ 2;

% Create direction sensor model --------------------------------------------------------------------
directionSensorModel = DirectionSensorModel();
directionSensorModel.sensorID = 2;

directionSensorModel.params.x.value = -1;
directionSensorModel.params.y.value = 0;
directionSensorModel.params.z.value = 0;

directionSensorModel.noiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);

% === Batch Algorithm ==============================================================================

% Set parameters to be estimated -------------------------------------------------------------------
projectileModel.params.CD.isEstimated = true;
earthModel.params.vWindx.isEstimated = true;

% Create batch estimator ---------------------------------------------------------------------------
estimator = BatchEstimator(projectileModelDynamics, integrator, { rangeSensorModel, directionSensorModel });

% Run batch algorithm ------------------------------------------------------------------------------
estimator.solve(estimateTime, sensorSampleRecord);


% ==================================================================================================
% === PLOTTING =====================================================================================
% ==================================================================================================

% Plot sensor measurements -------------------------------------------------------------------------

rangeSensorSampleRecord = sensorSampleRecord(rangeSensor.sensorID);

rangeSampleTimeHistory = rangeSensorSampleRecord.timeHistory;
rangeSampleHistory = rangeSensorSampleRecord.sampleHistory;

directionSensorSampleRecord = sensorSampleRecord(directionSensor.sensorID);

directionSampleTimeHistory = directionSensorSampleRecord.timeHistory;
azimuthSampleHistory = directionSensorSampleRecord.sampleHistory(1, :);
elevationSampleHistory = directionSensorSampleRecord.sampleHistory(2, :);

figure(1)

subplot(1, 3, 1)
plot(rangeSampleTimeHistory, rangeSampleHistory, 'kx')
xlabel("t (s)")
ylabel("R (m)")

subplot(1, 3, 2)
plot(directionSampleTimeHistory, rad2deg(azimuthSampleHistory), 'kx')
xlabel("t (s)")
ylabel("\theta (deg)")

subplot(1, 3, 3)
plot(directionSampleTimeHistory, rad2deg(elevationSampleHistory), 'kx')
xlabel("t (s)")
ylabel("\phi (deg)")

sgtitle("Sensor Measurements")

% Plot projectile trajectories ---------------------------------------------------------------------

nPlotPts = 250;
[~, truthStatePlotHistory] = truthPropagator.generateStateHistory(initSimTime, finalSimTime, nPlotPts);

figure(2)

plot3(truthStatePlotHistory(1, :), truthStatePlotHistory(2, :), -truthStatePlotHistory(3, :), 'k', "LineWidth", 1.5)
hold on
plot3(truthStatePlotHistory(1, end), truthStatePlotHistory(2, end), -truthStatePlotHistory(3, end), 'ko', "MarkerFaceColor", 'k', "HandleVisibility", "off")
for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    estimator.propagator.stateInterpolant = estimator.stateInterpolantList{ii};
    [~, nomStatePlotPts] = estimator.propagator.generateStateHistory(initSimTime, finalSimTime, nPlotPts);

    if ii < (Constants.MAX_ESTIMATOR_ITERATIONS + 1)
        plot3(nomStatePlotPts(1, :), nomStatePlotPts(2, :), -nomStatePlotPts(3, :), "Color", [1, 0.75, 0.75], "LineWidth", 0.5, "HandleVisibility", "off")
    else
        plot3(nomStatePlotPts(1, :), nomStatePlotPts(2, :), -nomStatePlotPts(3, :), 'r', "LineWidth", 1.5)
    end
end
hold off
grid on
xlabel("x (m)")
ylabel("y (m)")
zlabel("h (m)")
legend(["True", "Fit"])

sgtitle("Projectile Trajectory")

% Plot estimator convergence -----------------------------------------------------------------------

estimatedStateHistory = zeros(projectileModel.state.N_STATES, Constants.MAX_ESTIMATOR_ITERATIONS + 1);
for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    estimator.propagator.stateInterpolant = estimator.stateInterpolantList{ii};
    estimatedInitState = estimator.propagator.getState(estimateTime);

    estimatedStateHistory(:, ii) = estimatedInitState;
end

truthState = truthPropagator.getState(estimateTime);

figure(3)

subplot(1, 3, 1)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(4, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(4) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlabel("Iterations")
ylabel("v_x (m/s)")
legend(["Estimate", "Truth"])

subplot(1, 3, 2)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(5, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(5) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlabel("Iterations")
ylabel("v_y (m/s)")
legend(["Estimate", "Truth"])

subplot(1, 3, 3)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(6, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(6) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlabel("Iterations")
ylabel("v_z (m/s)")
legend(["Estimate", "Truth"])

sgtitle("Velocity Convergence")

% Plot measurment residuals ------------------------------------------------------------------------

priorRangeResidualRecord = estimator.measurementResidualRecordList{1}(rangeSensorModel.sensorID);

rangeResidualTimeHistory = priorRangeResidualRecord.timeHistory;
priorRangeResidualHistory = priorRangeResidualRecord.residualHistory;

postRangeResidualRecord = estimator.measurementResidualRecordList{end}(rangeSensorModel.sensorID);

postRangeResidualHistory = postRangeResidualRecord.residualHistory;

priorDirectionResidualRecord = estimator.measurementResidualRecordList{1}(directionSensorModel.sensorID);

directionResidualTimeHistory = priorDirectionResidualRecord.timeHistory;
priorAzimuthResidualHistory = priorDirectionResidualRecord.residualHistory(1, :);
priorElevationResidualHistory = priorDirectionResidualRecord.residualHistory(2, :);

postDirectionResidualRecord = estimator.measurementResidualRecordList{end}(directionSensorModel.sensorID);

postAzimuthResidualHistory = postDirectionResidualRecord.residualHistory(1, :);
postElevationResidualHistory = postDirectionResidualRecord.residualHistory(2, :);

figure(4)

subplot(1, 3, 1)
plot(rangeResidualTimeHistory, priorRangeResidualHistory, 'rx')
hold on
plot(rangeResidualTimeHistory, postRangeResidualHistory, 'kx')
hold off
xlabel("t (s)")
ylabel("\epsilon_R (m)")
legend(["Prefit", "Postfit"])

subplot(1, 3, 2)
plot(directionResidualTimeHistory, rad2deg(priorAzimuthResidualHistory), 'rx')
hold on
plot(directionResidualTimeHistory, rad2deg(postAzimuthResidualHistory), 'kx')
hold off
xlabel("t (s)")
ylabel("\epsilon_\theta (deg)")
legend(["Prefit", "Postfit"])

subplot(1, 3, 3)
plot(directionResidualTimeHistory, rad2deg(priorElevationResidualHistory), 'rx')
hold on
plot(directionResidualTimeHistory, rad2deg(postElevationResidualHistory), 'kx')
hold off
xlabel("t (s)")
ylabel("\epsilon_\phi (deg)")
legend(["Prefit", "Postfit"])

sgtitle("Measurement Residuals")

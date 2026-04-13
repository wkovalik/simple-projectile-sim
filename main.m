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
truthProjectile.state.vector = [0; 0; 39.5; 39.75];

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

rangeSensor.samplePeriod = 0.1;

rangeSensor.noiseCovar = 0.01 ^ 2;

% Create elevation sensor --------------------------------------------------------------------------
elevationSensor = ElevationSensor();
elevationSensor.sensorID = 2;

elevationSensor.params.x.value = -1;
elevationSensor.params.y.value = 0;

elevationSensor.samplePeriod = 0.01;

elevationSensor.noiseCovar = deg2rad(0.1) ^ 2;

% Create integrator --------------------------------------------------------------------------------
integrator = RK4Integrator();

integrator.stepPeriod = 0.01;

% Create truth propagator --------------------------------------------------------------------------
truthPropagator = Propagator(truthProjectileDynamics, integrator);

truthPropagator.addSensor(rangeSensor);
truthPropagator.addSensor(elevationSensor);

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
projectileModel.state.vector = [0; 0; 40; 40];
projectileModel.state.covar = diag([0.1; 0.1; 0.5; 0.5] .^ 2);

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

rangeSensorModel.noiseCovar = 0.01 ^ 2;

% Create elevation sensor model --------------------------------------------------------------------
elevationSensorModel = ElevationSensorModel();
elevationSensorModel.sensorID = 2;

elevationSensorModel.params.x.value = -1;
elevationSensorModel.params.y.value = 0;

elevationSensorModel.noiseCovar = deg2rad(0.1) ^ 2;

% === Batch Algorithm ==============================================================================

% Set parameters to be estimated -------------------------------------------------------------------
projectileModel.params.CD.isEstimated = true;
earthModel.params.vWindx.isEstimated = true;

% Create batch estimator ---------------------------------------------------------------------------
estimator = BatchEstimator(projectileModelDynamics, integrator, { rangeSensorModel, elevationSensorModel });

% Run batch algorithm ------------------------------------------------------------------------------
estimator.solve(estimateTime, sensorSampleRecord);

% TODO TOMORROW: 3-D with azimuth (directionSensor = elevation + azimuth). Storage. Kalman. Numerical Jacobian. Auto Jacobian. Comment.


% ==================================================================================================
% === PLOTTING =====================================================================================
% ==================================================================================================

% Plot sensor measurements -------------------------------------------------------------------------

rangeSensorSampleRecord = sensorSampleRecord(rangeSensor.sensorID);

rangeSampleTimeHistory = rangeSensorSampleRecord.timeHistory;
rangeSampleHistory = rangeSensorSampleRecord.sampleHistory;

elevationSensorSampleRecord = sensorSampleRecord(elevationSensor.sensorID);

elevationSampleTimeHistory = elevationSensorSampleRecord.timeHistory;
elevationSampleHistory = elevationSensorSampleRecord.sampleHistory;

figure(1)

subplot(1, 2, 1)
plot(rangeSampleTimeHistory, rangeSampleHistory, 'kx')
xlabel("t (s)")
ylabel("R (m)")
title("Range Measurements")

subplot(1, 2, 2)
plot(elevationSampleTimeHistory, rad2deg(elevationSampleHistory), 'kx')
xlabel("t (s)")
ylabel("\theta (deg)")
title("Elevation Measurements")

% Plot projectile trajectories ---------------------------------------------------------------------

nPlotPts = 250;
[~, truthStatePlotHistory] = truthPropagator.generateStateHistory(initSimTime, finalSimTime, nPlotPts);

figure(2)

plot(truthStatePlotHistory(1, :), truthStatePlotHistory(2, :), 'k', "LineWidth", 1.5)
hold on
for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    estimator.propagator.stateInterpolant = estimator.stateInterpolantList{ii};
    [~, nomStatePlotPts] = estimator.propagator.generateStateHistory(initSimTime, finalSimTime, nPlotPts);

    if ii < (Constants.MAX_ESTIMATOR_ITERATIONS + 1)
        plot(nomStatePlotPts(1, :), nomStatePlotPts(2, :), "Color", [1, 0.75, 0.75], "LineWidth", 0.5, "HandleVisibility", "off")
    else
        plot(nomStatePlotPts(1, :), nomStatePlotPts(2, :), 'r', "LineWidth", 1.5)
    end
end
hold off
xlabel("x (m)")
ylabel("y (m)")
legend(["True", "Fit"])
title("Projectile Trajectory")

% Plot estimator convergence -----------------------------------------------------------------------

estimatedStateHistory = zeros(projectileModel.state.N_STATES, Constants.MAX_ESTIMATOR_ITERATIONS + 1);
for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    estimator.propagator.stateInterpolant = estimator.stateInterpolantList{ii};
    estimatedInitState = estimator.propagator.getState(estimateTime);

    estimatedStateHistory(:, ii) = estimatedInitState;
end

truthState = truthPropagator.getState(estimateTime);

figure(3)

plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(3, :), 'x-', "LineWidth", 1.5)
hold on
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(4, :), 'x-', "LineWidth", 1.5)
hold off
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(3) * [1, 1], 'k', "LineWidth", 0.5)
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(4) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlabel("Iterations")
ylabel("v (m/s)")
legend(["v_x", "v_y"])
title("Fit State Convergence")

% Plot measurment residuals ------------------------------------------------------------------------

priorRangeResidualRecord = estimator.measurementResidualRecordList{1}(rangeSensorModel.sensorID);

priorRangeResidualTimeHistory = priorRangeResidualRecord.timeHistory;
priorRangeResidualHistory = priorRangeResidualRecord.residualHistory;

postRangeResidualRecord = estimator.measurementResidualRecordList{end}(rangeSensorModel.sensorID);

postRangeResidualTimeHistory = postRangeResidualRecord.timeHistory;
postRangeResidualHistory = postRangeResidualRecord.residualHistory;

priorElevationResidualRecord = estimator.measurementResidualRecordList{1}(elevationSensorModel.sensorID);

priorElevationResidualTimeHistory = priorElevationResidualRecord.timeHistory;
priorElevationResidualHistory = priorElevationResidualRecord.residualHistory;

postElevationResidualRecord = estimator.measurementResidualRecordList{end}(elevationSensorModel.sensorID);

postElevationResidualTimeHistory = postElevationResidualRecord.timeHistory;
postElevationResidualHistory = postElevationResidualRecord.residualHistory;

figure(4)

subplot(1, 2, 1)
plot(priorRangeResidualTimeHistory, priorRangeResidualHistory, 'rx')
hold on
plot(postRangeResidualTimeHistory, postRangeResidualHistory, 'kx')
hold off
xlabel("t (s)")
ylabel("\epsilon_R (m)")
legend(["Prefit", "Postfit"])
title("Range Residuals")

subplot(1, 2, 2)
plot(priorElevationResidualTimeHistory, priorElevationResidualHistory, 'rx')
hold on
plot(postElevationResidualTimeHistory, postElevationResidualHistory, 'kx')
hold off
xlabel("t (s)")
ylabel("\epsilon_\theta (deg)")
legend(["Prefit", "Postfit"])
title("Elevation Residuals")

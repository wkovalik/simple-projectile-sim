clear; clc; close all;

rng(0);

initSimTime = 0;   % (s) Initial time of truth simulation
finalSimTime = 7;  % (s) Final time of truth simulation

estimateTime = 0;  % (s) Time at which you want to estimate state


% ==================================================================================================
% === TRUTH TRAJECTORY =============================================================================
% ==================================================================================================

% Create truth projectile --------------------------------------------------------------------------
truthProjectile = Projectile();

% Define projectile params
truthProjectile.params.m.value  = 0.145;                 % (kg)
truthProjectile.params.S.value  = (pi / 4) * 0.075 ^ 2;  % (m^2)
truthProjectile.params.CD.value = 0.15;

% Define initial projectile state
truthProjectile.time = initSimTime;
truthProjectile.state.vector = [0; 0; 0; 39.5; 0; -39.75];  % (m; m; m; m/s; m/s; m/s)

% Create truth planet ------------------------------------------------------------------------------
truthEarth = Earth();

% Define planet parameters
truthEarth.params.vWindx.value = -5;  % (m/s)
truthEarth.params.vWindy.value = 10;  % (m/s)

% Create truth projectile dynamics -----------------------------------------------------------------
truthProjectileDynamics = ProjectileDynamics(truthProjectile, truthEarth);

% Create range sensor ------------------------------------------------------------------------------
rangeSensor = RangeSensor();
rangeSensor.sensorID = 1;  % Used to identify source of sensor measurements

% Define sensor parameters
rangeSensor.params.x.value = -1;  % (m)
rangeSensor.params.y.value = 0;   % (m)
rangeSensor.params.z.value = 0;   % (m)

% Define measurement sample period and covariance
rangeSensor.samplePeriod = 0.1;             % (s)
rangeSensor.noiseCovar = deg2rad(0.1) ^ 2;  % (m)^2

% Create direction sensor --------------------------------------------------------------------------
directionSensor = DirectionSensor();
directionSensor.sensorID = 2;  % Used to identify source of sensor measurements

% Define sensor parameters
directionSensor.params.x.value = -1;  % (m)
directionSensor.params.y.value = 0;   % (m)
directionSensor.params.z.value = 0;   % (m)

% Define measurement sample period and covariance
directionSensor.samplePeriod = 0.01;                          % (s)
directionSensor.noiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);  % (rad; rad)^2

% Create integrator --------------------------------------------------------------------------------
integrator = RK4Integrator();

integrator.stepPeriod = 0.01;  % (s)

% Create truth propagator --------------------------------------------------------------------------
truthPropagator = Propagator(truthProjectileDynamics, integrator);

% Add any sensors
truthPropagator.addSensor(rangeSensor);
truthPropagator.addSensor(directionSensor);

% Choose propagation option
truthPropagator.selectPropagator("stateWithSensors");  % Option to take sensor measurements during propagation

% Propagate truth trajectory -----------------------------------------------------------------------
truthPropagator.propagate(finalSimTime);

sensorSampleHistory = truthPropagator.sensorSampleHistory;  % Contains sensor measurement histories


% ==================================================================================================
% === ESTIMATION ===================================================================================
% ==================================================================================================

% Create projectile model --------------------------------------------------------------------------
projectileModel = Projectile();

% Define projectile params and covariance (only required for estimated parameters)
projectileModel.params.m.value  = 0.145;                 % (kg)
projectileModel.params.S.value  = (pi / 4) * 0.075 ^ 2;  % (m^2)
projectileModel.params.CD.value = 0.20;
projectileModel.params.CD.covar = 0.1 ^ 2;

% Define initial projectile state and covariance
projectileModel.time = initSimTime;
projectileModel.state.vector = [0; 0; 0; 40; 0; -40];                     % (m; m; m; m/s; m/s; m/s)
projectileModel.state.covar = diag([0.1; 0.1; 0.1; 0.5; 0.5; 0.5] .^ 2);  % (m; m; m; m/s; m/s; m/s)^2

% Create planet model ------------------------------------------------------------------------------
earthModel = Earth();

% Define planet parameters and covariance (only required for estimated parameters)
earthModel.params.vWindx.value = 0;       % (m/s)
earthModel.params.vWindy.value = 0;       % (m/s)
earthModel.params.vWindx.covar = 25 ^ 2;  % (m/s)^2
earthModel.params.vWindy.covar = 25 ^ 2;  % (m/s)^2

% Create projectile model dynamics -----------------------------------------------------------------
projectileModelDynamics = ProjectileDynamics(projectileModel, earthModel);

% Create range sensor model ------------------------------------------------------------------------
rangeSensorModel = RangeSensorModel();
rangeSensorModel.sensorID = 1;  % Must match ID of corresponding sensor

% Define sensor parameters
rangeSensorModel.params.x.value = -1;  % (m)
rangeSensorModel.params.y.value = 0;   % (m)
rangeSensorModel.params.z.value = 0;   % (m)

% Define measurement covariance
rangeSensorModel.noiseCovar = 0.01 ^ 2;  % (m)^2

% Create direction sensor model --------------------------------------------------------------------
directionSensorModel = DirectionSensorModel();
directionSensorModel.sensorID = 2;  % Must match ID of corresponding sensor

% Define sensor parameters
directionSensorModel.params.x.value = -1;  % (m)
directionSensorModel.params.y.value = 0;   % (m)
directionSensorModel.params.z.value = 0;   % (m)

% Define measurement covariance
directionSensorModel.noiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);  % (rad; rad)^2

% === Batch Algorithm ==============================================================================

% Set parameters to be estimated -------------------------------------------------------------------
projectileModel.params.CD.isEstimated = true;
earthModel.params.vWindx.isEstimated = true;
earthModel.params.vWindy.isEstimated = true;

% Create batch estimator ---------------------------------------------------------------------------
estimator = BatchEstimator(projectileModelDynamics, integrator, { rangeSensorModel, directionSensorModel });

% Run batch algorithm ------------------------------------------------------------------------------
estimator.solve(estimateTime, sensorSampleHistory);


% ==================================================================================================
% === PLOTTING =====================================================================================
% ==================================================================================================

% Plot sensor measurements -------------------------------------------------------------------------

% Get range sensor measurement history
rangeSensorSampleHistory = sensorSampleHistory(rangeSensor.sensorID);

rangeSampleTimeHistory = rangeSensorSampleHistory.timeHistory;
rangeSampleHistory = rangeSensorSampleHistory.sampleHistory;

% Get direction sensor measurement history
directionSensorSampleHistory = sensorSampleHistory(directionSensor.sensorID);

directionSampleTimeHistory = directionSensorSampleHistory.timeHistory;
azimuthSampleHistory = directionSensorSampleHistory.sampleHistory(1, :);
elevationSampleHistory = directionSensorSampleHistory.sampleHistory(2, :);

% Plot
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

% Get state history for truth trajectory
nPlotPts = 250;
[~, truthStatePlotHistory] = truthPropagator.generateStateHistory(initSimTime, finalSimTime, nPlotPts);

% Get state history for estimated trajectories (after each iteration of estimator)
estimatedStatePlotPts = cell(1, Constants.MAX_ESTIMATOR_ITERATIONS + 1);

for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    estimator.propagator.stateInterpolant = estimator.stateInterpolantHistory{ii};
    [~, estimatedStatePlotPts{ii}] = estimator.propagator.generateStateHistory(initSimTime, finalSimTime, nPlotPts);
end

% Plot
figure(2)

plot3(truthStatePlotHistory(1, :), truthStatePlotHistory(2, :), -truthStatePlotHistory(3, :), 'k', "LineWidth", 1.5)
hold on
plot3(truthStatePlotHistory(1, end), truthStatePlotHistory(2, end), -truthStatePlotHistory(3, end), 'ko', "MarkerFaceColor", 'k', "HandleVisibility", "off")
for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    if ii < (Constants.MAX_ESTIMATOR_ITERATIONS + 1)
        plot3(estimatedStatePlotPts{ii}(1, :), estimatedStatePlotPts{ii}(2, :), -estimatedStatePlotPts{ii}(3, :), "Color", [1, 0.75, 0.75], "LineWidth", 0.5, "HandleVisibility", "off")
    else
        plot3(estimatedStatePlotPts{ii}(1, :), estimatedStatePlotPts{ii}(2, :), -estimatedStatePlotPts{ii}(3, :), 'r', "LineWidth", 1.5)
    end
end
hold off
grid on
xlabel("x (m)")
ylabel("y (m)")
zlabel("h (m)")
legend(["True", "Fit"], "Location", "best")

sgtitle("Projectile Trajectory")

% Plot projectile velocity estimate convergence ----------------------------------------------------

% Get truth state
truthState = truthPropagator.getState(estimateTime);

% Get estimated states (after each iteration of estimator)
estimatedStateHistory = zeros(projectileModel.state.N_STATES, Constants.MAX_ESTIMATOR_ITERATIONS + 1);

for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    estimator.propagator.stateInterpolant = estimator.stateInterpolantHistory{ii};
    estimatedInitState = estimator.propagator.getState(estimateTime);

    estimatedStateHistory(:, ii) = estimatedInitState;
end

% Plot
figure(3)

subplot(1, 3, 1)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(4, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(4) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlim([0, Constants.MAX_ESTIMATOR_ITERATIONS])
xlabel("Iterations")
ylabel("v_x (m/s)")
legend(["Estimate", "Truth"])

subplot(1, 3, 2)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(5, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(5) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlim([0, Constants.MAX_ESTIMATOR_ITERATIONS])
xlabel("Iterations")
ylabel("v_y (m/s)")
legend(["Estimate", "Truth"])

subplot(1, 3, 3)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedStateHistory(6, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthState(6) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlim([0, Constants.MAX_ESTIMATOR_ITERATIONS])
xlabel("Iterations")
ylabel("v_z (m/s)")
legend(["Estimate", "Truth"])

sgtitle("Projectile Velocity Convergence")

% Plot parameter estimate convergence --------------------------------------------------------------

% Get estimated parameters (after each iteration of estimator)
estimatedCDHistory = zeros(1, Constants.MAX_ESTIMATOR_ITERATIONS + 1);
estimatedWindHistory = zeros(2, Constants.MAX_ESTIMATOR_ITERATIONS + 1);

for ii = 1:(Constants.MAX_ESTIMATOR_ITERATIONS + 1)
    estimatedCDHistory(ii) = estimator.estimatedParamHistory{ii}(1);
    estimatedWindHistory(:, ii) = estimator.estimatedParamHistory{ii}(2:3);
end

% Plot
figure(4)

subplot(1, 3, 1)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedCDHistory, 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthProjectile.params.CD.value * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlim([0, Constants.MAX_ESTIMATOR_ITERATIONS])
xlabel("Iterations")
ylabel("C_D")
legend(["Estimate", "Truth"])

subplot(1, 3, 2)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedWindHistory(1, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthEarth.params.vWindx.value * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlim([0, Constants.MAX_ESTIMATOR_ITERATIONS])
xlabel("Iterations")
ylabel("vWind_x (m/s)")
legend(["Estimate", "Truth"])

subplot(1, 3, 3)
plot(0:Constants.MAX_ESTIMATOR_ITERATIONS, estimatedWindHistory(2, :), 'rx-', "LineWidth", 1.5)
hold on
plot([0, Constants.MAX_ESTIMATOR_ITERATIONS], truthEarth.params.vWindy.value * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlim([0, Constants.MAX_ESTIMATOR_ITERATIONS])
xlabel("Iterations")
ylabel("vWind_y (m/s)")
legend(["Estimate", "Truth"])

sgtitle("Parameter Convergence")

% Plot measurment residuals ------------------------------------------------------------------------

% Get prefit range residuals
priorRangeResidualHistory = estimator.measurementResidualHistory{1}(rangeSensorModel.sensorID);

rangeResidualTimeHistory = priorRangeResidualHistory.timeHistory;
priorRangeResidualHistory = priorRangeResidualHistory.residualHistory;

% Get postfit range residuals
postRangeResidualHistory = estimator.measurementResidualHistory{end}(rangeSensorModel.sensorID);

postRangeResidualHistory = postRangeResidualHistory.residualHistory;

% Get prefit azimuth and elevation residuals
priorDirectionResidualHistory = estimator.measurementResidualHistory{1}(directionSensorModel.sensorID);

directionResidualTimeHistory = priorDirectionResidualHistory.timeHistory;
priorAzimuthResidualHistory = priorDirectionResidualHistory.residualHistory(1, :);
priorElevationResidualHistory = priorDirectionResidualHistory.residualHistory(2, :);

% Get postfit azimuth and elevation residuals
postDirectionResidualHistory = estimator.measurementResidualHistory{end}(directionSensorModel.sensorID);

postAzimuthResidualHistory = postDirectionResidualHistory.residualHistory(1, :);
postElevationResidualHistory = postDirectionResidualHistory.residualHistory(2, :);

% Plot
figure(5)

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

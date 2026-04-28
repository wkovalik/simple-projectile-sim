clear; clc; close all;

% TODO: This actually fails to converge within 10 iterations if CD0 = 0.20 +/- 0.10
% (It eventually converges within 15 iterations, but it still wanders everywhere in the beginning)
% Interestingly, it will converge if only 20 seconds of measurement history are included (versus
% the full 23-second history until ground impact). Investigate restarting estimator at apogee?

rng(0);

propagateTruthTrajectory();
runEstimator();
plotResults();

load("./log/trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
load("./log/sensorLog.mat", "measHistory");
load("./log/estimatorLog.mat", "output");



function propagateTruthTrajectory()
    fprintf("Propagating truth trajectory... ")

    % ----------------------------------------------------------------------------------------------
    
    % Create truth planet (default models)
    earth = Earth();
    
    % Set planet parameters
    earth.paramDefs.vWindx.value = -5;
    earth.paramDefs.vWindy.value = 10;
    
    earth.updateParams();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create truth projectile (default models)
    projectile = Projectile();
    
    % Set initial time and state
    projectile.stateDef.time = 0;
    projectile.stateDef.state = [0; 0; 0; 30; 0; -330];
    
    projectile.updateState();
    
    % Set projectile parameters
    projectile.paramDefs.CD.value = 0.15;
    
    projectile.updateParams();
    
    % Create projectile dynamics
    projectileDynamics = ProjectileDynamics(projectile, earth);
    
    % ----------------------------------------------------------------------------------------------
    
    % Create range finder sensor
    rangeSensor = RangeSensor();
    
    % Set sensor measurement properties
    rangeSensor.ID = 1;
    rangeSensor.samplePeriod = 0.05;
    rangeSensor.measNoiseCovar = 0.01 ^ 2;
    
    % Set sensor parameters
    rangeSensor.paramDefs.x.value = -1;
    
    rangeSensor.updateParams();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create direction finder sensor
    directionSensor = DirectionSensor();
    
    % Set sensor measurement properties
    directionSensor.ID = 2;
    directionSensor.samplePeriod = 0.05;
    directionSensor.measNoiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);
    
    % Set sensor parameters
    directionSensor.paramDefs.x.value = -1;
    
    directionSensor.updateParams();

    % ----------------------------------------------------------------------------------------------

    % Create accelerometer sensor
    accelerometerSensor = AccelerometerSensor(projectile, projectileDynamics);

    accelerometerSensor.ID = 3;
    accelerometerSensor.samplePeriod = 0.05;
    accelerometerSensor.measNoiseCovar = 0.01 ^ 2;
    
    % ----------------------------------------------------------------------------------------------
    
    % Create propagator (default integrator)
    propagator = Propagator(projectileDynamics);
    
    % Propagate truth trajectory (and take measurements along trajectory)
    propTime = 30;
    [trueTimeHistory, trueStateHistory, measHistory] = propagator.propagateWithSensors(propTime, { rangeSensor, directionSensor, accelerometerSensor });
    
    save("./log/trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
    save("./log/sensorLog.mat", "measHistory");

    fprintf("Done.\n\n")
end


function output = runEstimator()
    fprintf("Running estimator... \n")

    % ----------------------------------------------------------------------------------------------
    
    % Create planet model (default models)
    earthModel = Earth();
    
    % Set planet parameters and parameter covariances
    earthModel.paramDefs.vWindx.value = 0;
    earthModel.paramDefs.vWindx.covar = 50 ^ 2;
    earthModel.paramDefs.vWindx.isEstimated = true;
    
    earthModel.paramDefs.vWindy.value = 0;
    earthModel.paramDefs.vWindy.covar = 50 ^ 2;
    earthModel.paramDefs.vWindy.isEstimated = true;
    
    earthModel.updateParams();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create projectile model (default models)
    projectileModel = Projectile();
    
    % Set initial time, state, and state covariances
    projectileModel.time = 0;
    projectileModel.stateDef.state = [0; 0; 0; 30; 0; -330];
    projectileModel.stateDef.covar = diag([0.01; 0.01; 0.01; 0.5; 0.5; 5] .^ 2);  % TODO: Translate (V, az, el) with covars to (vx, vy, vz)
    
    projectileModel.updateState();
    
    % Set projectile parameters and parameter covariances
    projectileModel.paramDefs.CD.value = 0.165;
    projectileModel.paramDefs.CD.covar = 0.025 ^ 2;
    projectileModel.paramDefs.CD.isEstimated = true;
    
    projectileModel.updateParams();
    
    % Create projectile dynamics model
    projectileModelDynamics = ProjectileDynamics(projectileModel, earthModel);
    
    % ----------------------------------------------------------------------------------------------
    
    % Create range finder model
    rangeSensorModel = RangeSensor();
    
    % Set sensor measurement properties
    rangeSensorModel.ID = 1;
    rangeSensorModel.measNoiseCovar = 0.01 ^ 2;
    
    % Set sensor parameters
    rangeSensorModel.paramDefs.x.value = -1;
    
    rangeSensorModel.updateParams();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create direction finder model
    directionSensorModel = DirectionSensor();
    
    % Set sensor measurement properties
    directionSensorModel.ID = 2;
    directionSensorModel.measNoiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);
    
    % Set sensor parameters
    directionSensorModel.paramDefs.x.value = -1;
    
    directionSensorModel.updateParams();

    % ----------------------------------------------------------------------------------------------

    % Create accelerometer model
    accelerometerSensorModel = AccelerometerSensor(projectileModel, projectileModelDynamics);

    accelerometerSensorModel.ID = 3;
    accelerometerSensorModel.measNoiseCovar = 0.01 ^ 2;
    
    % ----------------------------------------------------------------------------------------------
    
    % Create estimator (default integrator)
    estimator = BatchEstimator(projectileModelDynamics, { rangeSensorModel, directionSensorModel, accelerometerSensorModel });
    
    % Run estimator
    load("./log/sensorLog.mat", "measHistory");

    % profile on -historysize 100000000
    output = estimator.solve(measHistory);
    % profile viewer
    % profile off
    
    save("./log/estimatorLog.mat", "output");

    fprintf("Done.\n\n")
end


function plotResults()
    fprintf("Plotting results... \n")

    load("./log/trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
    load("./log/sensorLog.mat", "measHistory");
    load("./log/estimatorLog.mat", "output");
    
    % ----------------------------------------------------------------------------------------------
    
    % Resample truth trajectory for plotting
    propTime = trueTimeHistory(end);
    plotTimeHistory = linspace(0, propTime, 250);
    plotTrueStateHistory = Utils.resampleStateHistory(trueTimeHistory, trueStateHistory, plotTimeHistory);
    
    % Resample nominal trajectories for plotting
    nIterations = length(output.iterationData) - 1;
    
    plotNomStateHistories = cell(1, nIterations + 1);
    for i = 1:(nIterations + 1)
        plotNomStateHistories{i} = Utils.resampleStateHistory(output.iterationData{i}.nomTimeHistory, output.iterationData{i}.nomStateHistory, plotTimeHistory);
    end
    
    % Get measurement histories
    rangeHistory = measHistory(2:end, measHistory(1, :) == 1);
    dirHistory = measHistory(2:end, measHistory(1, :) == 2);
    accHistory = measHistory(2:end, measHistory(1, :) == 3);
    
    % Get measurement residual histories
    priorMeasResiduals = output.iterationData{1}.measResidualHistory;
    priorRangeResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 1);
    priorDirResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 2);
    priorAccResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 3);
    
    postMeasResiduals = output.iterationData{end}.measResidualHistory;
    postRangeResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 1);
    postDirResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 2);
    postAccResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 3);
    
    % Get convergence histories
    CDHistory = output.paramHistory(1, :);
    CDStdDevHistory = output.paramCovarHistory(1, :) .^ 0.5;
    CDPlusHistory = CDHistory + CDStdDevHistory;
    CDMinusHistory = CDHistory - CDStdDevHistory;

    vWindxHistory = output.paramHistory(2, :);
    vWindxStdDevHistory = output.paramCovarHistory(5, :) .^ 0.5;
    vWindxPlusHistory = vWindxHistory + vWindxStdDevHistory;
    vWindxMinusHistory = vWindxHistory - vWindxStdDevHistory;

    vWindyHistory = output.paramHistory(3, :);
    vWindyStdDevHistory = output.paramCovarHistory(9, :) .^ 0.5;
    vWindyPlusHistory = vWindyHistory + vWindyStdDevHistory;
    vWindyMinusHistory = vWindyHistory - vWindyStdDevHistory;
    
    % ----------------------------------------------------------------------------------------------
    
    figure(1)
    
    plot3(plotTrueStateHistory(1, :), plotTrueStateHistory(2, :), -plotTrueStateHistory(3, :), 'k', "LineWidth", 1.5)
    hold on
    for i = 1:(nIterations + 1)
        plotNomStateHistory = plotNomStateHistories{i};
    
        if i == (nIterations + 1)
            plot3(plotNomStateHistory(1, :), plotNomStateHistory(2, :), -plotNomStateHistory(3, :), 'r', "LineWidth", 1.5)
        else
            plot3(plotNomStateHistory(1, :), plotNomStateHistory(2, :), -plotNomStateHistory(3, :), "Color", [1, 0.75, 0.75], "LineWidth", 0.5, "HandleVisibility", "off")
        end
    end
    hold off
    grid on
    xlabel("x (m)")
    ylabel("y (m)")
    zlabel("h (m)")
    legend(["True", "Fit"], "Location", "best")
    
    sgtitle("Projectile Trajectories")
    
    % ----------------------------------------------------------------------------------------------
    
    figure(2)
    
    subplot(2, 2, 1)
    plot(rangeHistory(1, :), rangeHistory(2, :), 'kx')
    xlabel("t (s)")
    ylabel("R (m)")
    
    subplot(2, 2, 2)
    plot(dirHistory(1, :), rad2deg(dirHistory(2, :)), 'kx')
    xlabel("t (s)")
    ylabel("\theta (deg)")
    
    subplot(2, 2, 3)
    plot(dirHistory(1, :), rad2deg(dirHistory(3, :)), 'kx')
    xlabel("t (s)")
    ylabel("\phi (deg)")

    subplot(2, 2, 4)
    plot(accHistory(1, :), accHistory(2, :), 'kx')
    xlabel("t (s)")
    ylabel("a_V (m/s^2)")
    
    sgtitle("Sensor Measurements")
    
    % ----------------------------------------------------------------------------------------------
    
    figure(3)
    
    subplot(2, 2, 1)
    plot(priorRangeResiduals(1, :), priorRangeResiduals(2, :), 'rx')
    hold on
    plot(postRangeResiduals(1, :), postRangeResiduals(2, :), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("R (m)")

    subplot(2, 2, 2)
    plot(priorDirResiduals(1, :), rad2deg(priorDirResiduals(2, :)), 'rx')
    hold on
    plot(postDirResiduals(1, :), rad2deg(postDirResiduals(2, :)), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("\theta (deg)")

    subplot(2, 2, 3)
    plot(priorDirResiduals(1, :), rad2deg(priorDirResiduals(3, :)), 'rx')
    hold on
    plot(postDirResiduals(1, :), rad2deg(postDirResiduals(3, :)), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("\phi (deg)")

    subplot(2, 2, 4)
    plot(priorAccResiduals(1, :), priorAccResiduals(2, :), 'rx')
    hold on
    plot(postAccResiduals(1, :), postAccResiduals(2, :), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("a_V (m/s^2)")
    
    legend(["Prefit", "Postfit"], "location", "best")
    
    sgtitle("Measurement Residuals")
    
    % ----------------------------------------------------------------------------------------------
    
    figure(4)
    
    subplot(1, 3, 1)
    patch([0:nIterations, flip(0:nIterations)], [CDPlusHistory, flip(CDMinusHistory)], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, CDHistory, 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("C_D")

    subplot(1, 3, 2)
    patch([0:nIterations, flip(0:nIterations)], [vWindxPlusHistory, flip(vWindxMinusHistory)], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindxHistory, 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_x (m/s)")

    subplot(1, 3, 3)
    patch([0:nIterations, flip(0:nIterations)], [vWindyPlusHistory, flip(vWindyMinusHistory)], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindyHistory, 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_y (m/s)")

    sgtitle("Convergence History")
end

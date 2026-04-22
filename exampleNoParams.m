clear; clc; close all;

rng(0);

propagateTruthTrajectory();
runEstimator();
plotResults();

load("trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
load("sensorLog.mat", "measHistory");
load("estimatorLog.mat", "output");



function propagateTruthTrajectory()
    fprintf("Propagating truth trajectory... ")

    % ----------------------------------------------------------------------------------------------
    
    % Create truth planet (default models and parameters)
    earth = Earth();
    
    % ----------------------------------------------------------------------------------------------

    % Create truth projectile (default models and parameters)
    projectile = Projectile();
    
    % Set initial time and state
    projectile.stateDef.time = 0;
    projectile.stateDef.state = [0; 0; 0; 30; 0; -330];
    
    projectile.updateState();
    
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

    % Create propagator (default integrator)
    propagator = Propagator(projectileDynamics);
    
    % Propagate truth trajectory (and take measurements along trajectory)
    propTime = 30;
    [trueTimeHistory, trueStateHistory, measHistory] = propagator.propagateWithSensors(propTime, { rangeSensor, directionSensor });
    
    save("trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
    save("sensorLog.mat", "measHistory");

    fprintf("Done.\n\n")
end


function output = runEstimator()
    fprintf("Running estimator... \n")

    % ----------------------------------------------------------------------------------------------

    % Create planet model (default models and parameters)
    earthModel = Earth();
    
    % ----------------------------------------------------------------------------------------------

    % Create projectile model (default models and parameters)
    projectileModel = Projectile();
    
    % Set initial time, state, and state covariances
    projectileModel.time = 0;
    projectileModel.stateDef.state = [0; 0; 0; 35; 5; -325];
    projectileModel.stateDef.covar = diag([0.01; 0.01; 0.01; 0.5; 0.5; 5] .^ 2);  % TODO: Translate (V, az, el) with covars to (vx, vy, vz)
    
    projectileModel.updateState();
    
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
    
    % Create estimator (default integrator)
    estimator = BatchEstimator(projectileModelDynamics, { rangeSensorModel, directionSensorModel });
    
    % Run estimator
    load("sensorLog.mat", "measHistory");

    % profile on -historysize 100000000
    output = estimator.solve(measHistory);
    % profile viewer
    % profile off
    
    save("estimatorLog.mat", "output");

    fprintf("Done.\n\n")
end


function plotResults()
    fprintf("Plotting results... \n")

    load("trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
    load("sensorLog.mat", "measHistory");
    load("estimatorLog.mat", "output");
    
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
    
    % Get measurement residual histories
    priorMeasResiduals = output.iterationData{1}.measResidualHistory;
    priorRangeResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 1);
    priorDirResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 2);
    
    postMeasResiduals = output.iterationData{end}.measResidualHistory;
    postRangeResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 1);
    postDirResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 2);
    
    % Get convergence histories
    vHistory = output.stateHistory(4:6, :);
    vStdDevHistory = output.stateCovarHistory([22, 29, 36], :) .^ 0.5;
    vPlusHistory = vHistory + vStdDevHistory;
    vMinusHistory = vHistory - vStdDevHistory;
    
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
    
    subplot(1, 3, 1)
    plot(rangeHistory(1, :), rangeHistory(2, :), 'kx')
    xlabel("t (s)")
    ylabel("R (m)")
    
    subplot(1, 3, 2)
    plot(dirHistory(1, :), rad2deg(dirHistory(2, :)), 'kx')
    xlabel("t (s)")
    ylabel("\theta (deg)")
    
    subplot(1, 3, 3)
    plot(dirHistory(1, :), rad2deg(dirHistory(3, :)), 'kx')
    xlabel("t (s)")
    ylabel("\phi (deg)")
    
    sgtitle("Sensor Measurements")
    
    % ----------------------------------------------------------------------------------------------
    
    figure(3)
    
    subplot(1, 3, 1)
    plot(priorRangeResiduals(1, :), priorRangeResiduals(2, :), 'rx')
    hold on
    plot(postRangeResiduals(1, :), postRangeResiduals(2, :), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("R (m)")
    
    subplot(1, 3, 2)
    plot(priorDirResiduals(1, :), rad2deg(priorDirResiduals(2, :)), 'rx')
    hold on
    plot(postDirResiduals(1, :), rad2deg(postDirResiduals(2, :)), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("\theta (deg)")
    
    subplot(1, 3, 3)
    plot(priorDirResiduals(1, :), rad2deg(priorDirResiduals(3, :)), 'rx')
    hold on
    plot(postDirResiduals(1, :), rad2deg(postDirResiduals(3, :)), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("\phi (deg)")
    
    legend(["Prefit", "Postfit"], "location", "best")
    
    sgtitle("Measurement Residuals")
    
    % ----------------------------------------------------------------------------------------------
    
    figure(4)
    
    subplot(1, 3, 1)
    patch([0:nIterations, flip(0:nIterations)], [vPlusHistory(1, :), flip(vMinusHistory(1, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vHistory(1, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("v_x (m/s)")

    subplot(1, 3, 2)
    patch([0:nIterations, flip(0:nIterations)], [vPlusHistory(2, :), flip(vMinusHistory(2, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vHistory(2, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("v_y (m/s)")

    subplot(1, 3, 3)
    patch([0:nIterations, flip(0:nIterations)], [vPlusHistory(3, :), flip(vMinusHistory(3, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vHistory(3, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("v_z (m/s)")

    sgtitle("Convergence History")
end

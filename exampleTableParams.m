clear; clc; close all;

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
    
    % Create truth planet
    earth = Earth();
    
    % Set planet models
    earth.atmosphereModel = "exponential";
    earth.windModel = "table";
    
    earth.updateModels();
    
    % Set planet parameters
    earth.paramDefs.vWindx.xValues = linspace(0, 1000, 6)';
    earth.paramDefs.vWindx.yValues = linspace(-5, 15, 6)';
    
    earth.paramDefs.vWindy.xValues = linspace(0, 1000, 6)';
    earth.paramDefs.vWindy.yValues = linspace(10, 0, 6)';
    
    earth.updateParams();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create truth projectile
    projectile = Projectile();
    
    % Set projectile time and state
    projectile.stateDef.time = 0;
    projectile.stateDef.state = [0; 0; 0; 30; 0; -330];
    
    projectile.updateState();
    
    % Set projectile models
    projectile.aeroModel = "table";
    
    projectile.updateModels();
    
    % Set projectile parameters
    projectile.paramDefs.CD.xValues = linspace(0, 1, 6)';
    projectile.paramDefs.CD.yValues = linspace(0.15, 0.15, 6)';
    
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
    
    % Create propagator (default integrator)
    propagator = Propagator(projectileDynamics);
    
    % Propagate truth trajectory (and take measurements along trajectory)
    propTime = 30;
    [trueTimeHistory, trueStateHistory, measHistory] = propagator.propagateWithSensors(propTime, { rangeSensor, directionSensor });
    
    save("./log/trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
    save("./log/sensorLog.mat", "measHistory");

    fprintf("Done.\n\n")
end


function output = runEstimator()
    fprintf("Running estimator... \n")

    % ----------------------------------------------------------------------------------------------

    % Create planet model
    earthModel = Earth();
    
    % Set planet models
    earthModel.atmosphereModel = "exponential";
    earthModel.windModel = "table";
    
    earthModel.updateModels();
    
    % Set planet parameters and parameter covariances
    earthModel.paramDefs.H.value = 9000;
    earthModel.paramDefs.H.covar = 1000 ^ 2;
    earthModel.paramDefs.H.isEstimated = true;
    
    earthModel.paramDefs.vWindx.xValues = [0; 500; 1000];
    earthModel.paramDefs.vWindx.yValues = [0; 0; 0];
    earthModel.paramDefs.vWindx.yCovars = [50; 50; 50] .^ 2;
    earthModel.paramDefs.vWindx.yIsEstimated = [true; true; true];
    
    earthModel.paramDefs.vWindy.xValues = [0; 500; 1000];
    earthModel.paramDefs.vWindy.yValues = [0; 0; 0];
    earthModel.paramDefs.vWindy.yCovars = [50; 50; 50] .^ 2;
    earthModel.paramDefs.vWindy.yIsEstimated = [true; true; true];
    
    earthModel.updateParams();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create projectile model
    projectileModel = Projectile();
    
    % Set initial time, state, and state covariances
    projectileModel.time = 0;
    projectileModel.stateDef.state = [0; 0; 0; 30; 0; -330];
    projectileModel.stateDef.covar = diag([0.01; 0.01; 0.01; 0.5; 0.5; 5] .^ 2);  % TODO: Translate (V, az, el) with covars to (vx, vy, vz)
    
    projectileModel.updateState();
    
    % Set projectile models
    projectileModel.aeroModel = "table";
    
    projectileModel.updateModels();
    
    % Set projectile parameters and parameter covariances
    projectileModel.paramDefs.CD.xValues = linspace(0, 1, 6)';
    projectileModel.paramDefs.CD.yValues = linspace(0.15, 0.15, 6)';
    projectileModel.paramDefs.CD.yIsEstimated = false(1, 6)';
    
    projectileModel.updateParams();
    
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
    
    % Get measurement residual histories
    priorMeasResiduals = output.iterationData{1}.measResidualHistory;
    priorRangeResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 1);
    priorDirResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 2);
    
    postMeasResiduals = output.iterationData{end}.measResidualHistory;
    postRangeResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 1);
    postDirResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 2);
    
    % Get convergence histories
    HHistory = output.paramHistory(1, :);
    HStdDevHistory = output.paramCovarHistory(1, :) .^ 0.5;
    HPlusHistory = HHistory + HStdDevHistory;
    HMinusHistory = HHistory - HStdDevHistory;
    
    vWindxHistory = output.paramHistory(2:4, :);
    vWindxStdDevHistory = output.paramCovarHistory([9, 17, 25], :) .^ 0.5;
    vWindxPlusHistory = vWindxHistory + vWindxStdDevHistory;
    vWindxMinusHistory = vWindxHistory - vWindxStdDevHistory;
    
    vWindyHistory = output.paramHistory(5:7, :);
    vWindyStdDevHistory = output.paramCovarHistory([33, 41, 49], :) .^ 0.5;
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
    
    patch([0:nIterations, flip(0:nIterations)], [HPlusHistory, flip(HMinusHistory)], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, HHistory, 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("H (m)")
    
    sgtitle("Convergence History")
    
    % ----------------------------------------------------------------------------------------------
    
    figure(5)
    
    subplot(3, 2, 1)
    patch([0:nIterations, flip(0:nIterations)], [vWindxPlusHistory(1, :), flip(vWindxMinusHistory(1, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindxHistory(1, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{x0} (m/s)")
    
    subplot(3, 2, 3)
    patch([0:nIterations, flip(0:nIterations)], [vWindxPlusHistory(2, :), flip(vWindxMinusHistory(2, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindxHistory(2, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{x1} (m/s)")
    
    subplot(3, 2, 5)
    patch([0:nIterations, flip(0:nIterations)], [vWindxPlusHistory(3, :), flip(vWindxMinusHistory(3, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindxHistory(3, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{x2} (m/s)")
    
    subplot(3, 2, 2)
    patch([0:nIterations, flip(0:nIterations)], [vWindyPlusHistory(1, :), flip(vWindyMinusHistory(1, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindyHistory(1, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{y0} (m/s)")
    
    subplot(3, 2, 4)
    patch([0:nIterations, flip(0:nIterations)], [vWindyPlusHistory(2, :), flip(vWindyMinusHistory(2, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindyHistory(2, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{y1} (m/s)")
    
    subplot(3, 2, 6)
    patch([0:nIterations, flip(0:nIterations)], [vWindyPlusHistory(3, :), flip(vWindyMinusHistory(3, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindyHistory(3, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{y2} (m/s)")
    
    sgtitle("Convergence History")
end

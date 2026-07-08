clear; clc; close all;

rng(0);

% TODO: Set maximum integrator step size!

propagateTruthTrajectory();
runEstimator("sequential");  % "sequential" or "sqrtsequential"
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

    earth.update();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create truth projectile
    projectile = Projectile();
    
    % Set projectile time and state
    projectile.stateDef.time = 0;
    projectile.stateDef.state = [0; 0; 0; 30; 0; -330; 120.0];
    
    % Set projectile models
    projectile.aeroModel = "table";
    
    % Set projectile parameters
    projectile.readPropsFromFile("data\csv\ANFinnerProps.csv");
    projectile.readAeroModelTablesFromFile("data\csv\ANFinnerAeroUniform.csv");
    
    projectile.update();
    
    % Create projectile dynamics
    projectileDynamics = ProjectileDynamics(projectile, earth);

    % ----------------------------------------------------------------------------------------------
    
    % Create roll gyro sensor
    rollGyroSensor = RollGyroSensor();

    % Set sensor measurement properties
    rollGyroSensor.ID = 1;
    rollGyroSensor.samplePeriod = 0.05;
    rollGyroSensor.measNoiseCovar = 0.0175 ^ 2;
    
    % ----------------------------------------------------------------------------------------------
    
    % Create propagator
    propagator = Propagator(projectileDynamics);
    propagator.integrator.stepPeriod = 0.01;
    
    % Propagate truth trajectory (and take measurements along trajectory)
    propTime = 30;
    [trueTimeHistory, trueStateHistory, measHistory] = propagator.propagateWithSensors(propTime, { rollGyroSensor });
    
    save("./log/trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
    save("./log/sensorLog.mat", "measHistory");

    fprintf("Done.\n\n")
end


function output = runEstimator(option)
    fprintf("Running %s estimator... \n", option)

    % ----------------------------------------------------------------------------------------------

    % Create planet model
    earthModel = Earth();
    
    % Set planet models
    earthModel.atmosphereModel = "exponential";
    
    % Set planet parameters and parameter covariances
    earthModel.paramDefs.H.value = 9000;
    earthModel.paramDefs.H.covar = 2500 ^ 2;
    earthModel.paramDefs.H.isEstimated = true;
    
    earthModel.update();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create projectile model
    projectileModel = Projectile();
    
    % Set initial time, state, and state covariances
    projectileModel.time = 0;
    projectileModel.stateDef.state = [0; 0; 0; 30; 0; -330; 120.0];
    % projectileModel.stateDef.covar = diag([0.01; 0.01; 0.01; 0.5; 0.5; 5; 6.2832] .^ 2);  % TODO: Translate (V, az, el) with covars to (vx, vy, vz)
    projectileModel.stateDef.covar = diag([0; 0; 0; 0.5; 0.5; 5; 6.2832] .^ 2);  % TODO: Translate (V, az, el) with covars to (vx, vy, vz)
    
    % Set projectile models
    projectileModel.aeroModel = "table";
    
    % Set projectile parameters
    projectileModel.readPropsFromFile("data\csv\ANFinnerProps.csv");
    projectileModel.readAeroModelTablesFromFile("data\csv\ANFinnerAeroUniform.csv");

    projectileModel.paramDefs.m.covar = 0.1 ^ 2;
    projectileModel.paramDefs.m.isConsidered = true;
    
    projectileModel.update();
    
    projectileModelDynamics = ProjectileDynamics(projectileModel, earthModel);
    
    % ----------------------------------------------------------------------------------------------
    
    % Create roll gyro sensor
    rollGyroSensorModel = RollGyroSensor();
    
    % Set sensor measurement properties
    rollGyroSensorModel.ID = 1;
    rollGyroSensorModel.measNoiseCovar = 0.0175 ^ 2;

    % ----------------------------------------------------------------------------------------------
    
    % Create estimator
    switch option
        case "sequential"
            estimator = SequentialConsiderEstimator(projectileModelDynamics, { rollGyroSensorModel });
        case "sqrtsequential"
            estimator = SqrtSequentialConsiderEstimator(projectileModelDynamics, { rollGyroSensorModel });
        otherwise
            error("Invalid estimator option.")
    end

    estimator.propagator.integrator.stepPeriod = 0.01;
    
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
    plotTimeHistory = linspace(0, propTime, 500);
    plotTrueStateHistory = Utils.resampleStateHistory(trueTimeHistory, trueStateHistory, plotTimeHistory);

    % Resample nominal trajectories for plotting
    nIterations = length(output.perIterationData) - 1;
    
    plotNomStateHistories = cell(1, nIterations + 1);
    for i = 1:(nIterations + 1)
        plotNomStateHistories{i} = Utils.resampleStateHistory(output.perIterationData{i}.nomTimeHistory, output.perIterationData{i}.nomStateHistory, plotTimeHistory);
    end

    % Get measurement histories
    rollGyroHistory = measHistory(2:end, measHistory(1, :) == 1);

    % Get measurement residual histories
    priorMeasResiduals = output.perIterationData{1}.measResidualHistory;
    priorGyroResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 1);

    postMeasResiduals = output.perIterationData{end}.measResidualHistory;
    postGyroResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 1);

    % Get convergence histories
    HIterations = output.iterations.params(1, :);
    HStdDevIterations = output.iterations.paramCovar(1, :) .^ 0.5;
    HPlusIterations = HIterations + HStdDevIterations;
    HMinusIterations = HIterations - HStdDevIterations;
    
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
    set(gca, "YDir", "reverse")
    xlabel("x (m)")
    ylabel("y (m)")
    zlabel("h (m)")
    legend(["True", "Fit"], "Location", "best")
    
    sgtitle("Projectile Trajectories")
    
    % ----------------------------------------------------------------------------------------------

    figure(2)

    subplot(2, 1, 1)
    plot(rollGyroHistory(1, :), rollGyroHistory(2, :) / (2 * pi), 'kx')
    xlabel("t (s)")
    ylabel("p (rev/s)")

    subplot(2, 1, 2)
    plot(priorGyroResiduals(1, :), priorGyroResiduals(2, :) / (2 * pi), 'rx')
    hold on
    plot(postGyroResiduals(1, :), postGyroResiduals(2, :) / (2 * pi), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("p (rev/s)")

    legend(["Prefit", "Postfit"], "location", "best")

    sgtitle("Measurements")

    % ----------------------------------------------------------------------------------------------
    
    figure(3)
    
    patch([0:nIterations, flip(0:nIterations)], [HPlusIterations, flip(HMinusIterations)], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, HIterations, 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("H (m)")
    
    sgtitle("Convergence History")
    

    % figure(2)
    % 
    % plot(plotTimeHistory, plotTrueStateHistory(7, :), 'k', "LineWidth", 1.5)
    % % plot(plotTrueStateHistory(1, :), -plotTrueStateHistory(3, :), 'k', "LineWidth", 1.5)
    % % plot(plotTimeHistory, rad2deg((plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5 ./ (plotTrueStateHistory(7, :) .^ 2 + plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5), 'k', "LineWidth", 1.5)
    % % plot(plotTimeHistory, (plotTrueStateHistory(7, :) .^ 2 + plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5 / 340.2824255476264, 'k', "LineWidth", 1.5)
    % grid on
end

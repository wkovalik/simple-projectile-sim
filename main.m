clear; clc; close all;

rng(0);

% TODO: Set maximum integrator step size!

propagateTruthTrajectory();
runEstimator("sqrtsequential");  % "batch", "squential", or "sqrtsequential"
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
    earthModel.paramDefs.H.covar = 1000 ^ 2;
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
    
    projectileModel.update();
    
    projectileModelDynamics = ProjectileDynamics(projectileModel, earthModel);
    
    % ----------------------------------------------------------------------------------------------
    
    % Create roll gyro sensor
    rollGyroSensorModel = RollGyroSensor();
    
    % Set sensor measurement properties
    rollGyroSensorModel.ID = 1;
    rollGyroSensorModel.measNoiseCovar = 0.0175 ^ 2;
    
    rollGyroSensorModel.update();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create estimator
    switch option
        case "batch"
            estimator = BatchEstimator(projectileModelDynamics, { rollGyroSensorModel });
        case "sequential"
            estimator = SequentialEstimator(projectileModelDynamics, { rollGyroSensorModel });
        case "sqrtsequential"
            estimator = SqrtSequentialEstimator(projectileModelDynamics, { rollGyroSensorModel });
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

    % ----------------------------------------------------------------------------------------------

    % Resample truth trajectory for plotting
    propTime = trueTimeHistory(end);
    plotTimeHistory = linspace(0, propTime, 500);
    plotTrueStateHistory = Utils.resampleStateHistory(trueTimeHistory, trueStateHistory, plotTimeHistory);

    % Get measurement histories
    rollGyroHistory = measHistory(2:end, measHistory(1, :) == 1);
    
    % ----------------------------------------------------------------------------------------------
    
    figure(1)
    
    plot3(plotTrueStateHistory(1, :), plotTrueStateHistory(2, :), -plotTrueStateHistory(3, :), 'k', "LineWidth", 1.5)
    grid on
    xlabel("x (m)")
    ylabel("y (m)")
    zlabel("h (m)")
    set(gca, "YDir", "reverse")
    
    % ----------------------------------------------------------------------------------------------

    figure(2)

    plot(rollGyroHistory(1, :), rollGyroHistory(2, :) / (2 * pi), 'kx')
    grid on
    xlabel("t (s)")
    ylabel("p (rev/s)")
    

    % figure(2)
    % 
    % plot(plotTimeHistory, plotTrueStateHistory(7, :), 'k', "LineWidth", 1.5)
    % % plot(plotTrueStateHistory(1, :), -plotTrueStateHistory(3, :), 'k', "LineWidth", 1.5)
    % % plot(plotTimeHistory, rad2deg((plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5 ./ (plotTrueStateHistory(7, :) .^ 2 + plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5), 'k', "LineWidth", 1.5)
    % % plot(plotTimeHistory, (plotTrueStateHistory(7, :) .^ 2 + plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5 / 340.2824255476264, 'k', "LineWidth", 1.5)
    % grid on
end

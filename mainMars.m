clear; clc; close all;

rng(0);

% TODO: Set maximum integrator step size!

% propagateTruthTrajectory();
% runEstimator("sqrtsequential");  % "sequential", or "sqrtsequential"
plotResults();

load("./log/trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
load("./log/sensorLog.mat", "measHistory");
load("./log/estimatorLog.mat", "output");



function propagateTruthTrajectory()
    fprintf("Propagating truth trajectory... ")

    % ----------------------------------------------------------------------------------------------
    
    % Create truth planet
    mars = Mars();

    % Set planet models
    mars.atmosphereModel = "exponential";
    mars.windModel = "table";
    mars.windModelKernel = "linear";

    % Set planet parameters
    mars.paramDefs.vWindx.xValues = linspace(0, 15000, 6)';
    mars.paramDefs.vWindx.yValues = linspace(0, 10, 6)';
    
    mars.paramDefs.vWindy.xValues = linspace(0, 15000, 6)';
    mars.paramDefs.vWindy.yValues = linspace(0, -15, 6)';

    mars.update();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create truth projectile
    projectile = Projectile();
    
    % Set projectile time and state
    projectile.stateDef.time = 0;
    projectile.stateDef.state = [0; 0; 0;
                                 0; 1.4801; 0;
                                 330.0; 0; 0;
                                 3.4558E+03; 0; 0];
    
    % Set projectile models
    projectile.aeroModel = "table";
    
    % Set projectile parameters
    projectile.readPropsFromFile("data\csv\ANFinnerProps.csv");
    projectile.readAeroModelTablesFromFile("data\csv\ANFinnerAeroUniform.csv");
    
    projectile.update();
    
    % Create projectile dynamics
    projectileDynamics = ProjectileDynamics(projectile, mars);

    % ----------------------------------------------------------------------------------------------
    
    % Create range finder sensor
    rangeSensor = RangeSensor();
    
    % Set sensor measurement properties
    rangeSensor.ID = 1;
    rangeSensor.samplePeriod = 0.05;
    rangeSensor.measNoiseCovar = 0.01 ^ 2;
    
    % ----------------------------------------------------------------------------------------------
    
    % Create direction finder sensor
    directionSensor = DirectionSensor();
    
    % Set sensor measurement properties
    directionSensor.ID = 2;
    directionSensor.samplePeriod = 0.05;
    directionSensor.measNoiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);

    % ----------------------------------------------------------------------------------------------
    
    % Create roll gyro sensor
    rollGyroSensor = RollGyroSensor();

    % Set sensor measurement properties
    rollGyroSensor.ID = 3;
    rollGyroSensor.samplePeriod = 0.05;
    rollGyroSensor.measNoiseCovar = 0.0175 ^ 2;
    
    % ----------------------------------------------------------------------------------------------
    
    % Create propagator
    propagator = Propagator(projectileDynamics);
    propagator.integrator.stepPeriod = 0.01;
    
    % Propagate truth trajectory (and take measurements along trajectory)
    propTime = 60;
    [trueTimeHistory, trueStateHistory, measHistory] = propagator.propagateWithSensors(propTime, { rangeSensor, directionSensor, rollGyroSensor });
    
    save("./log/trueTrajectoryLog.mat", "trueTimeHistory", "trueStateHistory");
    save("./log/sensorLog.mat", "measHistory");

    fprintf("Done.\n\n")
end


function output = runEstimator(option)
    fprintf("Running %s estimator... \n", option)

    % ----------------------------------------------------------------------------------------------

    % Create planet model
    marsModel = Mars();
    
    % Set planet models
    marsModel.atmosphereModel = "exponential";
    marsModel.windModel = "table";
    marsModel.windModelKernel = "linear";
    
    % Set planet parameters and parameter covariances
    marsModel.paramDefs.H.value = 12500;
    marsModel.paramDefs.H.covar = 2500 ^ 2;
    marsModel.paramDefs.H.isEstimated = true;

    marsModel.paramDefs.vWindx.xValues = [0; 7500; 15000];
    marsModel.paramDefs.vWindx.yValues = [0; 0; 0];
    marsModel.paramDefs.vWindx.yCovars = [50; 50; 50] .^ 2;
    marsModel.paramDefs.vWindx.yIsEstimated = [true; true; true];
    
    marsModel.paramDefs.vWindy.xValues = [0; 7500; 15000];
    marsModel.paramDefs.vWindy.yValues = [0; 0; 0];
    marsModel.paramDefs.vWindy.yCovars = [50; 50; 50] .^ 2;
    marsModel.paramDefs.vWindy.yIsEstimated = [true; true; true];
    
    marsModel.update();
    
    % ----------------------------------------------------------------------------------------------
    
    % Create projectile model
    projectileModel = Projectile();
    
    % Set initial time, state, and state covariances
    projectileModel.time = 0;
    projectileModel.stateDef.state = [0; 0; 0;
                                      0; 1.4801; 0;
                                      330.0; 0; 0;
                                      3.4558E+03; 0; 0];
    % projectileModel.stateDef.covar = diag([0.01; 0.01; 0.01; 0.5; 0.5; 5] .^ 2);  % TODO: Translate (V, az, el) with covars to (vx, vy, vz)
    projectileModel.stateDef.covar = diag([0; 0; 0; ...
                                           0; 0; 0; ...
                                           10; 0; 0; ...
                                           10.4720 * 10; 0.1745; 0.1745] .^ 2);  % TODO: Translate (V, az, el) with covars to (vx, vy, vz)
    
    % Set projectile models
    projectileModel.aeroModel = "table";
    
    % Set projectile parameters
    projectileModel.readPropsFromFile("data\csv\ANFinnerProps.csv");
    projectileModel.readAeroModelTablesFromFile("data\csv\ANFinnerAeroUniform.csv");
    
    projectileModel.update();
    
    projectileModelDynamics = ProjectileDynamics(projectileModel, marsModel);
    
    % ----------------------------------------------------------------------------------------------
    
    % Create range finder model
    rangeSensorModel = RangeSensor();
    
    % Set sensor measurement properties
    rangeSensorModel.ID = 1;
    rangeSensorModel.measNoiseCovar = 0.01 ^ 2;
    
    % ----------------------------------------------------------------------------------------------
    
    % Create direction finder model
    directionSensorModel = DirectionSensor();
    
    % Set sensor measurement properties
    directionSensorModel.ID = 2;
    directionSensorModel.measNoiseCovar = diag(deg2rad([0.1; 0.1]) .^ 2);
    
    % ----------------------------------------------------------------------------------------------
    
    % Create roll gyro sensor
    rollGyroSensorModel = RollGyroSensor();
    
    % Set sensor measurement properties
    rollGyroSensorModel.ID = 3;
    rollGyroSensorModel.measNoiseCovar = 0.0175 ^ 2;
    
    % ----------------------------------------------------------------------------------------------
    
    % Create estimator
    switch option
        case "sequential"
            estimator = SequentialEstimator(projectileModelDynamics, { rangeSensorModel, directionSensorModel, rollGyroSensorModel });
        case "sqrtsequential"
            estimator = SqrtSequentialEstimator(projectileModelDynamics, { rangeSensorModel, directionSensorModel, rollGyroSensorModel });
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
    rangeHistory = measHistory(2:end, measHistory(1, :) == 1);
    dirHistory = measHistory(2:end, measHistory(1, :) == 2);
    gyroHistory = measHistory(2:end, measHistory(1, :) == 3);

    % Get measurement residual histories
    priorMeasResiduals = output.perIterationData{1}.measResidualHistory;
    priorRangeResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 1);
    priorDirResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 2);
    priorGyroResiduals = priorMeasResiduals(2:end, priorMeasResiduals(1, :) == 3);

    postMeasResiduals = output.perIterationData{end}.measResidualHistory;
    postRangeResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 1);
    postDirResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 2);
    postGyroResiduals = postMeasResiduals(2:end, postMeasResiduals(1, :) == 3);

    % Get convergence histories
    HIterations = output.iterations.params(1, :);
    HStdDevIterations = output.iterations.paramCovar(1, :) .^ 0.5;
    HPlusIterations = HIterations + HStdDevIterations;
    HMinusIterations = HIterations - HStdDevIterations;

    vWindxIterations = output.iterations.params(2:4, :);
    vWindxStdDevIterations = output.iterations.paramCovar([9, 17, 25], :) .^ 0.5;
    vWindxPlusIterations = vWindxIterations + vWindxStdDevIterations;
    vWindxMinusIterations = vWindxIterations - vWindxStdDevIterations;
    
    vWindyIterations = output.iterations.params(5:7, :);
    vWindyStdDevIterations = output.iterations.paramCovar([33, 41, 49], :) .^ 0.5;
    vWindyPlusIterations = vWindyIterations + vWindyStdDevIterations;
    vWindyMinusIterations = vWindyIterations - vWindyStdDevIterations;
    
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
    plot(gyroHistory(1, :), gyroHistory(2, :) / (2 * pi), 'kx')
    xlabel("t (s)")
    ylabel("p (rev/s)")

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
    plot(priorGyroResiduals(1, :), priorGyroResiduals(2, :) / (2 * pi), 'rx')
    hold on
    plot(postGyroResiduals(1, :), postGyroResiduals(2, :) / (2 * pi), 'kx')
    hold off
    xlabel("t (s)")
    ylabel("p (rev/s)")

    legend(["Prefit", "Postfit"], "location", "best")

    sgtitle("Measurement Residuals")
    
    % ----------------------------------------------------------------------------------------------
    
    figure(4)
    
    patch([0:nIterations, flip(0:nIterations)], [HPlusIterations, flip(HMinusIterations)], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, HIterations, 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("H (m)")
    
    sgtitle("Convergence History")

    % ----------------------------------------------------------------------------------------------
    
    figure(5)
    
    subplot(3, 2, 1)
    patch([0:nIterations, flip(0:nIterations)], [vWindxPlusIterations(1, :), flip(vWindxMinusIterations(1, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindxIterations(1, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{x0} (m/s)")
    
    subplot(3, 2, 3)
    patch([0:nIterations, flip(0:nIterations)], [vWindxPlusIterations(2, :), flip(vWindxMinusIterations(2, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindxIterations(2, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{x1} (m/s)")
    
    subplot(3, 2, 5)
    patch([0:nIterations, flip(0:nIterations)], [vWindxPlusIterations(3, :), flip(vWindxMinusIterations(3, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindxIterations(3, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{x2} (m/s)")
    
    subplot(3, 2, 2)
    patch([0:nIterations, flip(0:nIterations)], [vWindyPlusIterations(1, :), flip(vWindyMinusIterations(1, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindyIterations(1, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{y0} (m/s)")
    
    subplot(3, 2, 4)
    patch([0:nIterations, flip(0:nIterations)], [vWindyPlusIterations(2, :), flip(vWindyMinusIterations(2, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindyIterations(2, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{y1} (m/s)")
    
    subplot(3, 2, 6)
    patch([0:nIterations, flip(0:nIterations)], [vWindyPlusIterations(3, :), flip(vWindyMinusIterations(3, :))], [1, 0.8, 0.8], "EdgeColor", "none")
    hold on
    plot(0:nIterations, vWindyIterations(3, :), 'rx-', "LineWidth", 1.5)
    hold off
    box on
    xlim([1, nIterations])
    xlabel("Iterations")
    ylabel("vWind_{y2} (m/s)")
    
    sgtitle("Convergence History")

    % figure(2)
    % 
    % plot(plotTimeHistory, plotTrueStateHistory(7, :), 'k', "LineWidth", 1.5)
    % % plot(plotTrueStateHistory(1, :), -plotTrueStateHistory(3, :), 'k', "LineWidth", 1.5)
    % % plot(plotTimeHistory, rad2deg((plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5 ./ (plotTrueStateHistory(7, :) .^ 2 + plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5), 'k', "LineWidth", 1.5)
    % % plot(plotTimeHistory, (plotTrueStateHistory(7, :) .^ 2 + plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5 / 340.2824255476264, 'k', "LineWidth", 1.5)
    % grid on
end

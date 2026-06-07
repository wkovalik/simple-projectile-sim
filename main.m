clear; clc; close all;


rng(0);

% ----------------------------------------------------------------------------------------------
    
% Create truth planet (default models)
earth = Earth();

earth.update();

% ----------------------------------------------------------------------------------------------

% Create truth projectile (default models)
projectile = Projectile();

% Set initial time and state
projectile.stateDef.time = 0;
projectile.stateDef.state = [0; 0; -3.048; 0; 0.017953; 0; 1029.0; 0; 0; 120.0; 1.0; 0];

% Set projectile models
projectile.aeroModel = "table";

% Set projectile parameters
projectile.readPropsFromFile("data\csv\ANFinnerProps.csv");
projectile.readAeroModelTablesFromFile("data\csv\ANFinnerAeroUniform.csv");

projectile.update();

% Create projectile dynamics
projectileDynamics = ProjectileDynamics(projectile, earth);

% ----------------------------------------------------------------------------------------------

% Create propagator
propagator = Propagator(projectileDynamics);
propagator.integrator.stepPeriod = 0.001;

% Propagate truth trajectory (and take measurements along trajectory)
propTime = 4;
[trueTimeHistory, trueStateHistory] = propagator.propagate(propTime);

fprintf("Done.\n\n")

% ----------------------------------------------------------------------------------------------

% Resample truth trajectory for plotting
propTime = trueTimeHistory(end);
plotTimeHistory = linspace(0, propTime, 250);
plotTrueStateHistory = Utils.resampleStateHistory(trueTimeHistory, trueStateHistory, plotTimeHistory);

% ----------------------------------------------------------------------------------------------

figure(1)

plot3(plotTrueStateHistory(1, :), plotTrueStateHistory(2, :), -plotTrueStateHistory(3, :), 'k', "LineWidth", 1.5)
grid on
xlabel("x (m)")
ylabel("y (m)")
zlabel("h (m)")
set(gca, "YDir", "reverse")

% % plot(plotTimeHistory, (plotTrueStateHistory(7, :) .^ 2 + plotTrueStateHistory(8, :) .^ 2 + plotTrueStateHistory(9, :) .^ 2) .^ 0.5 / 340.2824255476264, 'k', "LineWidth", 1.5)
% % plot(plotTimeHistory, plotTrueStateHistory(12, :), 'k', "LineWidth", 1.5)
% plot(plotTrueStateHistory(2, :), -plotTrueStateHistory(3, :), 'k', "LineWidth", 1.5)
% grid on

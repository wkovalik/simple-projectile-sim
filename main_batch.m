clear; clc; close all;

rng(0);

initSimTime = 0;
finalSimTime = 7;

% ==================================================================================================
% === TRUE TRAJECTORY ==============================================================================
% ==================================================================================================

trueState0 = [0; 0; 39.5; 39.75];

trueProjectile.m = 0.145;
trueProjectile.S = (pi / 4) * 0.075 ^ 2;
trueProjectile.CD = 0.15;
trueEarth.g = 9.81;
trueEarth.rho = 1.225;
trueEarth.vWindx = -5;

[timePts, trueStatePts] = ode45( ...
    @(t, state) projectileModel(t, state, trueProjectile, trueEarth), ...
    [initSimTime, finalSimTime], ...
    trueState0 ...
);

trueStateTraj = griddedInterpolant(timePts, trueStatePts, "linear");

% ==================================================================================================
% === MEASUREMENTS =================================================================================
% ==================================================================================================

% === Range Measurements ===========================================================================

trueRangeSensor.samplePd = 0.1;
trueRangeSensor.x = -1;
trueRangeSensor.y = 0;
trueRangeSensor.covar = 0.01 ^ 2;

rangeSampleTimes = initSimTime:trueRangeSensor.samplePd:finalSimTime;

nRangeSamples = length(rangeSampleTimes);
rangeSamples = zeros(1, nRangeSamples);

for i = 1:length(rangeSampleTimes)
    sampleTime = rangeSampleTimes(i);
    rangeSample = rangeSensorSample(sampleTime, trueStateTraj(sampleTime), trueRangeSensor);

    rangeSamples(i) = rangeSample;
end

% === Elevation Measurements =======================================================================

trueElevSensor.samplePd = 0.01;
trueElevSensor.x = -1;
trueElevSensor.y = 0;
trueElevSensor.covar = deg2rad(0.1) ^ 2;

elevSampleTimes = initSimTime:trueElevSensor.samplePd:finalSimTime;

nElevSamples = length(elevSampleTimes);
elevSamples = zeros(1, nElevSamples);

for i = 1:length(elevSampleTimes)
    sampleTime = elevSampleTimes(i);
    elevSample = elevSensorSample(sampleTime, trueStateTraj(sampleTime), trueElevSensor);

    elevSamples(i) = elevSample;
end

% === Combined Measurement History =================================================================

% Sensor ID enum: 1 = range sensor, 2 = elevation sensor

rangeSampleRecord = [ 1 * ones(1, length(rangeSampleTimes));
                      rangeSampleTimes;
                      rangeSamples ];
elevSampleRecord = [ 2 * ones(1, length(elevSampleTimes));
                     elevSampleTimes;
                     elevSamples ];

sampleRecord = [rangeSampleRecord, elevSampleRecord];

sampleRecord = sortrows(sampleRecord', [2, 1])';  % Sort based on sample epoch, then sensor ID

nSamples = size(sampleRecord, 2);
finalSampleTime = sampleRecord(2, end);

% ==================================================================================================
% === ESTIMATOR ====================================================================================
% ==================================================================================================

nStates = 4;
nEstStates = 4;
nEstParams = 2;

estStateIndices = 1:4;
estParamsIndices = (1:nEstParams) + nStates;
estIndices = [estStateIndices, estParamsIndices];


nomProjectile.m = 0.145;
nomProjectile.S = (pi / 4) * 0.075 ^ 2;
nomProjectile.CD = 0.20;
nomEarth.g = 9.81;
nomEarth.rho = 1.225;
nomEarth.vWindx = 0;

nomRangeSensor.x = -1;
nomRangeSensor.y = 0;
nomRangeSensor.covar = 0.01 ^ 2;
nomRangeSensor.invCovar = 1 / nomRangeSensor.covar;

nomElevSensor.x = -1;
nomElevSensor.y = 0;
nomElevSensor.covar = deg2rad(0.1) ^ 2;
nomElevSensor.invCovar = 1 / nomElevSensor.covar;

% TODO: Rename 0 to init

priorState0 = [0; 0; 40; 40];
priorEstParams = [nomProjectile.CD; nomEarth.vWindx];

priorEstStateP0 = diag([0.1; 0.1; 0.5; 0.5] .^ 2);
priorEstParamsP = diag([0.1; 25] .^ 2);

priorP0 = blkdiag(priorEstStateP0, priorEstParamsP);


% === Batch Algorithm ==============================================================================

% TODO: In actual implementation, would sort measurement samples here into separate range sample
% matrix, elevatio sample matrix, etc. Still would feed full sampleMatrix into estimator since it is
% time sorted, though. Corresponding residual matrices will be populated in estimator loops

maxIters = 5;

nomStateTrajs = cell(1, maxIters + 1);

measResiduals = cell(1, 2);  % 2 = num sensors
measResiduals{1} = zeros(maxIters + 1, nRangeSamples);
measResiduals{2} = zeros(maxIters + 1, nElevSamples);


nomState0 = priorState0;
nomEstParams = priorEstParams;

priorDev0 = zeros(nEstStates + nEstParams, 1);

for ii = 1:maxIters

    % === Nominal Trajectory =======================================================================

    stateSTM0 = eye(nStates);
    stateSTM0 = stateSTM0(:);

    paramSTM0 = zeros(nStates, nEstParams);
    paramSTM0 = paramSTM0(:);

    augNomState0 = [nomState0; stateSTM0; paramSTM0];

    [timePts, augNomStatePts] = ode45( ...
        @(t, augState) projectileModelWithSTM(t, augState, nomProjectile, nomEarth), ...
        [initSimTime, finalSampleTime], ...
        augNomState0 ...
    );

    nomStatePts = augNomStatePts(:, 1:nStates);
    stateSTMPts = augNomStatePts(:, (nStates + 1):(nStates ^ 2 + nStates));
    paramSTMPts = augNomStatePts(:, (nStates ^ 2 + nStates + 1):end);

    nomStateTraj = griddedInterpolant(timePts, nomStatePts, "linear");
    nomStateTrajs{ii} = nomStateTraj;

    stateSTMTraj = griddedInterpolant(timePts, stateSTMPts, "linear");
    paramSTMTraj = griddedInterpolant(timePts, paramSTMPts, "linear");

    % === Batch Update =============================================================================

    Lambda0 = inv(priorP0);
    N0 = priorP0 \ priorDev0;

    sampleCounts = ones(1, 2);

    for i = 1:nSamples
        sensorID = sampleRecord(1, i);
        sampleTime = sampleRecord(2, i);
        sampleMeas = sampleRecord(3, i);

        nomState = nomStateTraj(sampleTime)';

        stateSTM = stateSTMTraj(sampleTime)';
        stateSTM = reshape(stateSTM, [nStates, nStates]);

        paramSTM = paramSTMTraj(sampleTime)';
        paramSTM = reshape(paramSTM, [nStates, nEstParams]);

        switch sensorID
            case 1
                nomMeas = rangeMeasModel(sampleTime, nomState, nomRangeSensor);
                RInv = nomRangeSensor.invCovar;

                stateH = rangeMeasStateJacobian(sampleTime, nomState, nomRangeSensor);

                paramH = zeros(1, 2);
                paramH(1, 1) = rangeMeasDragJacobian(sampleTime, nomState, nomRangeSensor);
                paramH(1, 2) = rangeMeasWindJacobian(sampleTime, nomState, nomRangeSensor);
                
            case 2
                nomMeas = elevMeasModel(sampleTime, nomState, nomElevSensor);
                RInv = nomElevSensor.invCovar;

                stateH = elevMeasStateJacobian(sampleTime, nomState, nomElevSensor);
                
                paramH = zeros(1, 2);
                paramH(1, 1) = elevMeasDragJacobian(sampleTime, nomState, nomRangeSensor);
                paramH(1, 2) = elevMeasWindJacobian(sampleTime, nomState, nomRangeSensor);
        end

        measResidual = sampleMeas - nomMeas;
        
        measResiduals{sensorID}(ii, sampleCounts(sensorID)) = measResidual;
        sampleCounts(sensorID) = sampleCounts(sensorID) + 1;

        stateHtilde = stateH * stateSTM;
        paramHtilde = stateH * paramSTM + paramH;

        Htilde = [stateHtilde, paramHtilde];

        nextLambda0 = Htilde' * RInv * Htilde;
        nextN0 = Htilde' * RInv * measResidual;

        Lambda0 = Lambda0 + nextLambda0;
        N0 = N0 + nextN0;
    end
    
    % Compute post-measurement deviations
    postDev0 = Lambda0 \ N0;
    
    postStateDev0 = postDev0(1:nEstStates);
    postParamsDev = postDev0((nEstStates + 1):end);
    
    % Update nominal initial states
    nomState0 = nomState0 + postStateDev0;
    
    % Update nominal parameters
    nomEstParams = nomEstParams + postParamsDev;

    nomProjectile.CD = nomEstParams(1);
    nomEarth.vWindx = nomEstParams(2);

    % Update pre-measurement deviations
    priorDev0 = priorDev0 - postDev0;

end

% === Postfit Trajectory ===========================================================================

[timePts, nomStatePts] = ode45( ...
    @(t, state) projectileModel(t, state, nomProjectile, nomEarth), ...
    [initSimTime, finalSimTime], ...
    nomState0 ...
);

nomStateTraj = griddedInterpolant(timePts, nomStatePts, "linear");
nomStateTrajs{end} = nomStateTraj;

% === Postfit Residuals ============================================================================

sampleCounts = ones(1, 2);

for i = 1:nSamples
    sensorID = sampleRecord(1, i);
    sampleTime = sampleRecord(2, i);
    sampleMeas = sampleRecord(3, i);

    nomState = nomStateTraj(sampleTime)';

    switch sensorID
        case 1
            nomMeas = rangeMeasModel(sampleTime, nomState, nomRangeSensor);
        case 2
            nomMeas = elevMeasModel(sampleTime, nomState, nomElevSensor);
    end

    measResidual = sampleMeas - nomMeas;
    
    measResiduals{sensorID}(end, sampleCounts(sensorID)) = measResidual;
    sampleCounts(sensorID) = sampleCounts(sensorID) + 1;
end


% ==================================================================================================
% === PLOTTING =====================================================================================
% ==================================================================================================

% === Range Measurements ===========================================================================

figure(1)

subplot(1, 2, 1)
plot(rangeSampleTimes, rangeSamples, 'kx')
xlabel("t (s)")
ylabel("R (m)")
title("Range Measurements")

subplot(1, 2, 2)
plot(elevSampleTimes, rad2deg(elevSamples), 'kx')
xlabel("t (s)")
ylabel("\theta (deg)")
title("Elevation Measurements")

% === Projectile Trajectories ======================================================================

nPlotPts = 250;
tPlotPts = linspace(initSimTime, finalSimTime, nPlotPts);
trueStatePlotPts = trueStateTraj(tPlotPts)';

figure(2)
plot(trueStatePlotPts(1, :), trueStatePlotPts(2, :), 'k', "LineWidth", 1.5)
hold on
for ii = 1:maxIters
    nomStateTraj = nomStateTrajs{ii};
    nomStatePlotPts = nomStateTraj(tPlotPts)';

    if ii < maxIters
        plot(nomStatePlotPts(1, :), nomStatePlotPts(2, :), 'r', "LineWidth", 0.5, "HandleVisibility", "off")
    else
        plot(nomStatePlotPts(1, :), nomStatePlotPts(2, :), 'r', "LineWidth", 1.5)
    end
end
hold off
xlabel("x (m)")
ylabel("y (m)")
legend(["True", "Fit"])
title("Projectile Trajectory")

% === Estimate Convergence =========================================================================

nomState0Pts = zeros(4, maxIters);
for ii = 1:(maxIters + 1)
    nomStateTraj = nomStateTrajs{ii};
    nomState0 = nomStateTraj(0)';
    nomState0Pts(:, ii) = nomState0;
end

figure(3)
plot(0:maxIters, nomState0Pts(3, :), 'x-', "LineWidth", 1.5)
hold on
plot(0:maxIters, nomState0Pts(4, :), 'x-', "LineWidth", 1.5)
hold off
hold on
plot([0, maxIters], trueState0(3) * [1, 1], 'k', "LineWidth", 0.5)
plot([0, maxIters], trueState0(4) * [1, 1], 'k', "LineWidth", 0.5)
hold off
xlabel("Iterations")
ylabel("v (m/s)")
legend(["v_x", "v_y"])
title("Fit State Convergence")

% === Postfit Residuals ============================================================================

rangeResiduals = measResiduals{1};
elevResiduals = measResiduals{2};

figure(4)

subplot(1, 2, 1)
plot(rangeSampleTimes, rangeResiduals(1, :), 'rx')
hold on
plot(rangeSampleTimes, rangeResiduals(end, :), 'kx')
hold off
xlabel("t (s)")
ylabel("\epsilon_R (m)")
legend(["Prefit", "Postfit"])
title("Range Residuals")

subplot(1, 2, 2)
plot(elevSampleTimes, rad2deg(elevResiduals(1, :)), 'rx')
hold on
plot(elevSampleTimes, rad2deg(elevResiduals(end, :)), 'kx')
hold off
xlabel("t (s)")
ylabel("\epsilon_\theta (deg)")
legend(["Prefit", "Postfit"])
title("Elevation Residuals")


state0Error = trueState0 - nomState0;
fprintf("%.9f\n", state0Error(1))
fprintf("%.9f\n", state0Error(2))
fprintf("%.9f\n", state0Error(3))
fprintf("%.9f\n", state0Error(4))


% ==================================================================================================
% === FUNCTIONS ====================================================================================
% ==================================================================================================

function [Fx, Fy] = projectileGravModel(~, ~, projectile, earth)

m = projectile.m;
g = earth.g;

Fx = 0;
Fy = -m * g;

end


function [Fx, Fy] = projectileAeroModel(~, state, projectile, earth)

vx = state(3);
vy = state(4);

S = projectile.S;
CD = projectile.CD;
rho = earth.rho;
vWindx = earth.vWindx;

vInfx = vx - vWindx;
VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;

Fx = -(rho * S * CD / 2) * VInf * vInfx;
Fy = -(rho * S * CD / 2) * VInf * vy;

end


function stateDeriv = projectileModel(t, state, projectile, earth)

vx = state(3);
vy = state(4);

m = projectile.m;

[FGravx, FGravy] = projectileGravModel(t, state, projectile, earth);
[FAerox, FAeroy] = projectileAeroModel(t, state, projectile, earth);

Fx = FGravx + FAerox;
Fy = FGravy + FAeroy;

xDeriv = vx;
yDeriv = vy;
vxDeriv = Fx / m;
vyDeriv = Fy / m;

stateDeriv = [xDeriv; yDeriv; vxDeriv; vyDeriv];

end


function A = projectileStateJacobian(~, state, projectile, earth)

vx = state(3);
vy = state(4);

m = projectile.m;
S = projectile.S;
CD = projectile.CD;
rho = earth.rho;
vWindx = earth.vWindx;

vInfx = vx - vWindx;
VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;

A = zeros(4, 4);

A(1, 3) = 1;
A(3, 3) = -(rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
A(4, 3) = -(rho * S * CD / (2 * m)) * (vInfx * vy / VInf);

A(2, 4) = 1;
A(3, 4) = -(rho * S * CD / (2 * m)) * (vInfx * vy / VInf);
A(4, 4) = -(rho * S * CD / (2 * m)) * (vy ^ 2 / VInf + VInf);

end


function A = projectileDragJacobian(~, state, projectile, earth)

vx = state(3);
vy = state(4);

m = projectile.m;
S = projectile.S;
rho = earth.rho;
vWindx = earth.vWindx;

vInfx = vx - vWindx;
VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;

A = zeros(4, 1);

A(3, 1) = -(rho * S / (2 * m)) * VInf * vInfx;
A(4, 1) = -(rho * S / (2 * m)) * VInf * vy;

end


function A = projectileWindJacobian(~, state, projectile, earth)

vx = state(3);
vy = state(4);

m = projectile.m;
S = projectile.S;
CD = projectile.CD;
rho = earth.rho;
vWindx = earth.vWindx;

vInfx = vx - vWindx;
VInf = (vInfx ^ 2 + vy ^ 2) ^ 0.5;

A = zeros(4, 1);

A(3, 1) = (rho * S * CD / (2 * m)) * (vInfx ^ 2 / VInf + VInf);
A(4, 1) = (rho * S * CD / (2 * m)) * (vInfx * vy / VInf);

end


function augStateDeriv = projectileModelWithSTM(t, augState, projectile, earth)

state = augState(1:4);

stateDeriv = projectileModel(t, state, projectile, earth);


stateSTM = augState(5:20);
stateSTM = reshape(stateSTM, [4, 4]);

stateA = projectileStateJacobian(t, state, projectile, earth);

stateSTMDeriv = stateA * stateSTM;
stateSTMDeriv = stateSTMDeriv(:);


paramSTM = augState(21:end);
paramSTM = reshape(paramSTM, [4, 2]);

paramA = zeros(4, 2);
paramA(:, 1) = projectileDragJacobian(t, state, projectile, earth);
paramA(:, 2) = projectileWindJacobian(t, state, projectile, earth);

paramSTMDeriv = stateA * paramSTM + paramA;
paramSTMDeriv = paramSTMDeriv(:);


augStateDeriv = [stateDeriv; stateSTMDeriv; paramSTMDeriv];

end


function R = rangeSensorSample(~, state, sensor)

x = state(1);
y = state(2);

xStation = sensor.x;
yStation = sensor.y;
sigmaStation = sensor.covar ^ 0.5;

dxStation = x - xStation;
dyStation = y - yStation;

trueR = (dxStation ^ 2 + dyStation ^ 2) ^ 0.5;

R = trueR + sigmaStation * randn();

end


function R = rangeMeasModel(~, state, sensor)

x = state(1);
y = state(2);

xStation = sensor.x;
yStation = sensor.y;

dxStation = x - xStation;
dyStation = y - yStation;

R = (dxStation ^ 2 + dyStation ^ 2) ^ 0.5;

end


function H = rangeMeasStateJacobian(~, state, sensor)

x = state(1);
y = state(2);

xStation = sensor.x;
yStation = sensor.y;

dxStation = x - xStation;
dyStation = y - yStation;

R = (dxStation ^ 2 + dyStation ^ 2) ^ 0.5;

H = zeros(1, 4);

H(1, 1) = dxStation / R;

H(1, 2) = dyStation / R;

end


function H = rangeMeasDragJacobian(~, ~, ~)

H = 0;

end


function H = rangeMeasWindJacobian(~, ~, ~)

H = 0;

end


function theta = elevSensorSample(~, state, sensor)

x = state(1);
y = state(2);

xStation = sensor.x;
yStation = sensor.y;
sigmaStation = sensor.covar ^ 0.5;

dxStation = x - xStation;
dyStation = y - yStation;

trueTheta = atan2(dyStation, dxStation);

theta = trueTheta + sigmaStation * randn();

end


function theta = elevMeasModel(~, state, sensor)

x = state(1);
y = state(2);

xStation = sensor.x;
yStation = sensor.y;

dxStation = x - xStation;
dyStation = y - yStation;

theta = atan2(dyStation, dxStation);

end


function H = elevMeasStateJacobian(~, state, sensor)

x = state(1);
y = state(2);

xStation = sensor.x;
yStation = sensor.y;

dxStation = x - xStation;
dyStation = y - yStation;

denom = 1 + (dyStation / dxStation) ^ 2;

H = zeros(1, 4);

H(1, 1) = -(dyStation / dxStation ^ 2) / denom;

H(1, 2) = (1 / dxStation) / denom;

end


function H = elevMeasDragJacobian(~, ~, ~)

H = 0;

end


function H = elevMeasWindJacobian(~, ~, ~)

H = 0;

end

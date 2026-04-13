clear; clc; close all;


% Select parameters to estimate
params.P1.isEstimated = true;
params.P2.isEstimated = false;
params.P3.isEstimated = true;


% Build list of estimated parameter names (act as keys for Jacobian function dictionary)
estimatedParams = {};

paramNames = fieldnames(params);
for i = 1:length(paramNames)
    paramName = paramNames{i};

    if params.(paramName).isEstimated
        estimatedParams(end + 1) = { paramName };  %#ok<SAGROW> Suppress preallocation warning
    end
end

estimatedParams = string(estimatedParams);


% Use Jacobian function dictionary to build Jacobian for estimated parameters
jacobianMap = dictionary( ...
    "P1", @computeA1, ...
    "P2", @computeA2, ...
    "P3", @computeA3 ...
);

N_PROJECTILE_STATES = 3;
N_ESTIMATED_PARAMS = length(estimatedParams);

A = zeros(N_PROJECTILE_STATES, N_ESTIMATED_PARAMS);

for i = 1:N_ESTIMATED_PARAMS
    estimatedParam = estimatedParams(i);

    A(:, i) = feval(jacobianMap(estimatedParam));
end


function A = computeA1()
    A = [1; 10; 100];
end

function A = computeA2()
    A = [2; 20; 200];
end

function A = computeA3()
    A = [3; 30; 300];
end

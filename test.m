clear; clc; close all;

projectile.m = Param();
projectile.m.value = 1;
projectile.m.covar = 0.1;
projectile.m.isEstimated = true;

projectile.S = Param();
projectile.S.value = 1;
projectile.S.covar = 0;
projectile.S.isEstimated = false;

projectile.CD = Param();
projectile.CD.value = 0.1;
projectile.CD.covar = 0.01;
projectile.CD.isEstimated = true;

% projectile.paramNames = { projectile.m.name, ...
%                           projectile.S.name, ...
%                           projectile.CD.name };

paramNames = fieldnames(projectile);

estimatedParams = {};
for i = 1:length(paramNames)
    paramName = paramNames{i};

    if projectile.(paramName).isEstimated
        estimatedParams = [estimatedParams, projectile.(paramName)];
    end
end

estimatedParams
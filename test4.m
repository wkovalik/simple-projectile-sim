clear; clc; close all;


A1.sensorID = 1;
A1.sampleTime = 0.1;
A1.sampledMeasurement = 1.001;

A2.sensorID = 1;
A2.sampleTime = 0.15;
A2.sampledMeasurement = 2.001;

A3.sensorID = 1;
A3.sampleTime = 0.2;
A3.sampledMeasurement = 3.001;

A4.sensorID = 1;
A4.sampleTime = 0.25;
A4.sampledMeasurement = 4.001;

A5.sensorID = 2;
A5.sampleTime = 0.1;
A5.sampledMeasurement = [1.001; 1.001; 1.001];

A6.sensorID = 2;
A6.sampleTime = 0.25;
A6.sampledMeasurement = [2.001; 2.001; 2.001];

A = { A1, A2, A3, A4, A5, A6 };

times = zeros(1, length(A));
sensorIDs = zeros(1, length(A));
for i = 1:length(A)
    times(i) = A{i}.sampleTime;
    sensorIDs(i) = A{i}.sensorID;
end

[~, sortIndices] = sortrows([times; sensorIDs]', [1, 2]);

sortedA = cell(1, length(A));
for i = 1:length(A)
    sortedA{i} = A{sortIndices(i)};
end
clear; clc; close all;


sensor1.sensorID = 1;
sensor2.sensorID = 2;

sensor1.N_MEASUREMENTS = 1;
sensor2.N_MEASUREMENTS = 3;

sensorList = { sensor1, sensor2 };


sensorSampleRecord = initializeSensorSampleRecord(sensorList);


sensor1.sampleTime = 0.1;
sensor1.sampledMeasurement = 1.001;

sensorSampleRecord = addSample(sensorSampleRecord, sensor1);


sensor2.sampleTime = 0.1;
sensor2.sampledMeasurement = [1.001; 1.001; 1.001];

sensorSampleRecord = addSample(sensorSampleRecord, sensor2);


sensor1.sampleTime = 0.15;
sensor1.sampledMeasurement = 2.001;

sensorSampleRecord = addSample(sensorSampleRecord, sensor1);


sensor1.sampleTime = 0.2;
sensor1.sampledMeasurement = 3.001;

sensorSampleRecord = addSample(sensorSampleRecord, sensor1);


sensor1.sampleTime = 0.25;
sensor1.sampledMeasurement = 4.001;

sensorSampleRecord = addSample(sensorSampleRecord, sensor1);


sensor2.sampleTime = 0.25;
sensor2.sampledMeasurement = [2.001; 2.001; 2.001];

sensorSampleRecord = addSample(sensorSampleRecord, sensor2);


for i = 1:length(sensorList)
    sensor = sensorList{i};

    nSamples = sensorSampleRecord(sensor.sensorID).nSamples;
    sensorSampleRecord(sensor.sensorID).timeHistory((nSamples + 1):end) = [];
    sensorSampleRecord(sensor.sensorID).sampleHistory(:, (nSamples + 1):end) = [];
end


% sensorSampleRecord(sensor1.sensorID)
% sensorSampleRecord(sensor2.sensorID)

finalSampleTime = getFinalSampleTime(sensorList, sensorSampleRecord);

chronSensorSampleRecord = createChronSensorSampleRecord(sensorSampleRecord);



function sensorSampleRecord = initializeSensorSampleRecord(sensorList)
    sensorSampleRecord = dictionary();

    for i = 1:length(sensorList)
        sensor = sensorList{i};

        initRecord.nSamples = 0;
        initRecord.timeHistory = zeros(1, Constants.HISTORY_BUFFER_LEN);
        initRecord.sampleHistory = zeros(sensor.N_MEASUREMENTS, Constants.HISTORY_BUFFER_LEN);

        sensorSampleRecord = insert(sensorSampleRecord, sensor.sensorID, initRecord);
    end
end

function sensorSampleRecord = addSample(sensorSampleRecord, sensor)
    nSamples = sensorSampleRecord(sensor.sensorID).nSamples;
    nSamples = nSamples + 1;
    
    sensorSampleRecord(sensor.sensorID).nSamples = nSamples;
    sensorSampleRecord(sensor.sensorID).timeHistory(nSamples) = sensor.sampleTime;
    sensorSampleRecord(sensor.sensorID).sampleHistory(:, nSamples) = sensor.sampledMeasurement;
end

function chronSensorSampleRecord = createChronSensorSampleRecord(sensorSampleRecord)
    sensorIDs = keys(sensorSampleRecord);

    nTotalSamples = 0;
    for i = 1:length(sensorIDs)
        sensorID = sensorIDs(i);
        nSamples = sensorSampleRecord(sensorID).nSamples;

        nTotalSamples = nTotalSamples + nSamples;
    end

    concatSensorSampleRecord = cell(1, nTotalSamples);
    sampleCount = 0;
    for i = 1:length(sensorIDs)
        sensorID = sensorIDs(i);
        nSamples = sensorSampleRecord(sensorID).nSamples;

        for j = 1:nSamples
            sampleCount = sampleCount + 1;
            
            concatSensorSampleRecord{sampleCount}.sensorID = sensorID;
            concatSensorSampleRecord{sampleCount}.sampleTime = sensorSampleRecord(sensorID).timeHistory(j);
            concatSensorSampleRecord{sampleCount}.sampledMeasurement = sensorSampleRecord(sensorID).sampleHistory(:, j);
        end
    end

    sampleTimes = zeros(1, nTotalSamples);
    sensorIDs = zeros(1, nTotalSamples);
    for i = 1:nTotalSamples
        sampleTimes(i) = concatSensorSampleRecord{i}.sampleTime;
        sensorIDs(i) = concatSensorSampleRecord{i}.sensorID;
    end

    [~, sortIndices] = sortrows([sampleTimes; sensorIDs]', [1, 2]);

    chronSensorSampleRecord = cell(1, nTotalSamples);
    for i = 1:nTotalSamples
        chronSensorSampleRecord{i} = concatSensorSampleRecord{sortIndices(i)};
    end
end

function finalSampleTime = getFinalSampleTime(sensorList, sensorSampleRecord)
    finalSampleTime = 0;
    for i = 1:length(sensorList)
        sensor = sensorList{i};
        sensorFinalSampleTime = sensorSampleRecord(sensor.sensorID).timeHistory(end);
    
        if sensorFinalSampleTime > finalSampleTime
            finalSampleTime = sensorFinalSampleTime;
        end
    end
end

classdef Utils
    methods (Static)
        function stateHistory2 = resampleStateHistory(timeHistory1, stateHistory1, timeHistory2)
            stateInterpolant1 = griddedInterpolant( ...
                timeHistory1,                       ...
                stateHistory1',                     ...
                Constants.DEFAULT_INTERP_METHOD,    ...
                Constants.DEFAULT_EXTRAP_METHOD     ...
            );

            stateHistory2 = stateInterpolant1(timeHistory2)';
        end


        function nMeasMax = findLongestSensorMeas(sensorArray)
            nMeasMax = 0;

            for i = 1:length(sensorArray)
                sensor = sensorArray{i};
                
                if sensor.nMeas > nMeasMax
                    nMeasMax = sensor.nMeas;
                end
            end
        end


        function isTooCoarse = isIntegratorStepTooCoarse(sensorArray, integratorStepPeriod)
            isTooCoarse = false;

            for i = 1:length(sensorArray)
                sensor = sensorArray{i};

                if mod(sensor.samplePeriod, integratorStepPeriod) ~= 0
                    isTooCoarse = true;
                    break
                end
            end
        end
    end
end
classdef Validator
    methods (Static)
        function validValue = validateType(value, validType)
            % Checks if a value is of type validType

            if ~isa(validType, "string")
                error("Specified type must be a string.")
            end

            if isa(value, validType)
                validValue = value;
            else
                error("Invalid input. Must be type %s.", validType)
            end
        end

        function validStruct = validateStructArrayTypes(struct, validType)
            % Checks if all elements of a struct are of type validType (hence the name struct
            % "array", since all elements must be of a single type)

            fieldNames = fieldnames(struct);
            
            for i = 1:length(fieldNames)
                Validator.validateType(struct.(fieldNames{i}), validType);
            end

            validStruct = struct;
        end

        function validArray = validateSize(array, validSize)
            % Checks if an array has an (m x n) size validSize

            validSize = Validator.validateType(validSize, "double");
            
            if isequal(size(array), validSize)
                validArray = array;
            else
                error("Invalid input size. Must be size (%0.f, %0.f).", validSize(1), validSize(2))
            end
        end
    end
end
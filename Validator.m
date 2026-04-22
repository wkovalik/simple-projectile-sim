classdef Validator
    methods (Static)
        function validValue = validateType(value, validTypes)
            % Checks if a value is of any type in validTypes

            if ~isa(validTypes, "string")
                error("Specified type(s) must be a string.")
            end
            
            if isscalar(validTypes)
                if isa(value, validTypes)
                    validValue = value;
                else
                    error("Invalid input. Must be type %s.", validTypes)
                end

            else
                for i = 1:length(validTypes)
                    if isa(value, validTypes(i))
                        validValue = value;
                        break
                    end
                    
                    if i == length(validTypes)
                        errorStr = "Invalid input type. Must be one of the following options: [";
                        errorStr = strcat(errorStr, sprintf("'%s', ", validTypes{1:(end - 1)}));
                        errorStr = strcat(errorStr, sprintf("'%s'].", validTypes{end}));
        
                        error(errorStr)
                    end
                end
            end
        end


        function validStruct = validateFieldTypes(struct, validTypes)
            % Checks if all fields of a struct are of any type in validTypes

            fieldNames = fieldnames(struct);
            
            for i = 1:length(fieldNames)
                Validator.validateType(struct.(fieldNames{i}), validTypes);
            end

            validStruct = struct;
        end


        function validArray = validateLength(array, validLen)
            % Checks if an array has a length validLen (i.e., a (validLen, 1) size)

            Validator.validateType(validLen, "double");
            
            if length(array) == validLen
                validArray = array;
            else
                error("Invalid input length. Must be length %0.f.", validLen)
            end
        end


        function validArray = validateSize(array, validSize)
            % Checks if an array has an (m x n) size validSize

            Validator.validateType(validSize, "double");
            
            if isequal(size(array), validSize)
                validArray = array;
            else
                errorStr = "Invalid input size. Must be size (";
                if length(validSize) > 1
                    errorStr = strcat(errorStr, sprintf("%0.f, ", validSize(1:(end - 1))));
                end
                errorStr = strcat(errorStr, sprintf("%0.f).", validSize(end)));

                error(errorStr)
            end
        end

        
        function validStr = validateString(str, validStrs)
            % Checks if a string is inside an array of valid strings

            Validator.validateType(str, "string");
            Validator.validateType(validStrs, "string");

            if any(strcmp(str, validStrs))
                validStr = str;
            else
                errorStr = "Invalid input string. Must be one of the following options: [";
                if length(validStrs) > 1
                    errorStr = strcat(errorStr, sprintf("'%s', ", validStrs{1:(end - 1)}));
                end
                errorStr = strcat(errorStr, sprintf("'%s'].", validStrs{end}));

                error(errorStr)
            end
        end
    end
end
classdef FilesCycler
    %FILESCYCLER receives the parent directory of a file structure and it
    %opens the files ona at the time to analyze them. Saving results,
    %cnfigurations, and others will be implemented in future versions of
    %the class or in children classes.
    %   Detailed explanation goes here
    
    properties
        Directory char
    end
    properties (SetAccess = 'private')
        FileStructure
    end
    
    methods
        function obj = untitled3(inputArg1,inputArg2)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


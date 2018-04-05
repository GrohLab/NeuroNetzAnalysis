classdef MouseExperiment < handle
    %MOUSEEXPERIMENT Parent class dealing with the organization and
    %processing (?) of the experiments
    %   Detailed explanation goes here
    
    properties (SetAccess = 'private')
        % Unchangable properties of the experiment
        AnimalName
        ExperimentDate
        AnimalBirthday
        Experimenter
        Gender
    end
    
    properties (Dependent)
        PSTH
    end
    
    
    methods (Abstract)
        function obj = MouseExperiment(inputArg1,inputArg2)
            %MOUSEEXPERIMENT
            %   The constructor of this class imports the information from
            %   different types of file, from the workspace, or from manual
            %   input.
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


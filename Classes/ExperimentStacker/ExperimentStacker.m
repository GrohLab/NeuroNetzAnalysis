classdef ExperimentStacker < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'protected')
        Conditions
        Triggers
        DataDir char
        SamplingFrequency
        KSexe KSrunner
    end
    
    methods
        function obj = ExperimentStacker(DataDir)
            %EXPERIMENTSTACKER Construct an instance of this class using
            %the directory where all the experiment relevant files are
            %stored and either loads the data or creates the appropiated
            %variables and saves them.
            %   obj = ExperimentStacker(DataDir)
            if ~exist(DataDir,'dir')
                fprintf(1,'The given folder doesn''t exist.\n(%s)\n',DataDir)
                str =...
                    input('Open a GUI to select the correct one? (y/n): ',...
                    's');
                if ~isItYes(str)
                    obj.delete;
                    return
                end
                DataDir = uigetdir(pwd,...
                    'Please select the experiment folder');
            end
            obj.DataDir = DataDir;
        end
        
        function obj = prepareFilesForKS(obj)
            
        end
        
        function disp(obj)
            fprintf(1,'Experiment folder:\n%s',obj.DataDir)
        end
    end
end



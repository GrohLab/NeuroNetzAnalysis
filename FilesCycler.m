classdef FilesCycler
    %FILESCYCLER receives the parent directory of a file structure and it
    %opens the files ona at the time to analyze them. Saving results,
    %cnfigurations, and others will be implemented in future versions of
    %the class or in children classes.
    %   Detailed explanation goes here
    
    properties
        Directory char
        Function function_handler
    end
    properties (SetAccess = 'private')
        FileStructure
    end
    
    methods
        function obj = FilesCycler(ParentDir,FunctionHandler)
            %FILESCYCLER Construct an instance of this class
            %   Detailed explanation goes here
            if exist('ParentDir','var') && exist(ParentDir,'dir')
                obj.Directory = ParentDir;
            else
                fprintf('Invalid directory. Object not created!\n')
                delete(obj)
                return
            end
            if exist('FunctionHandler','var') && isa(FunctionHandler,'function_handler')
                obj.Function = FunctionHandler;
            else
                fprintf('Invalid function handler. Try setting it up correctly\n')
            end
        end
        
        function iok = validateFiles(obj,fileIdentifier)
            fileStruct = obj.FileStructure;
            
        end
        
        function fileStruct = get.FileStructure(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            currentContents = dir(obj.Directory);
            for cf = 1:
            
        end
    end
    
    methods (Static,Access = 'private')
        function recursiveAccess(directory,fileID)
            
        end
    end
end


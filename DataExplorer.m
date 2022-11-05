classdef DataExplorer < handle
    %{
    DATAEXPLORER This class handles the experiment files and searches for
    produces simple trigger outputs such as peri-stimulus time histograms,
    raster plots, triggered average of continuous signals, etc.
    
    This is the first version constructed from prototypes.
    
    Emilio Isaias-Camacho @ GrohLab 2019
    %}   
    
    properties (Access = 'private')
        configStruct struct
    end
    
    methods
        function obj = DataExplorer(inputArg1,inputArg2)
           %UNTITLED Construct an instance of this class
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


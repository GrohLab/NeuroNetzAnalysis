classdef Stacks
    %STACKS This class will be usefull for creating the PSTH's for which
    %ever 'discrete' signal (e.g. spikes, triggers, step functions),
    %triggered average signals from the 'continuous' waveforms, or raster
    %plots from the multi unit recordings. The main function of this class
    %is to give a support for the flexibility required to 'explore' the
    %data using such basic analyses. 
    %   Instance creation, methods and properties attributes. Inputs:
    %   signals from the workspace, file(s), interface with Kilosort are
    %   some of the possibilities for initiating an instance.
    %
    %   The methods will be rather simple and even abstract as the PSTH,
    %   raster, and triggered average are planned as 'independent' or child
    %   objects from this class.
    %
    %   The properties attributes should be given a round of thoughts
    %   because they might make the object (or instance) vulnerable to
    %   unintentional and error-creating changes. Finally, it is worth
    %   mentioning that programatically, it is, at least for me, easier to
    %   solve some technical issues if the class is a child from the MATLAB
    %   handle class. The object changes the values of itself without
    %   needing an auxiliary instance.
    properties
        Property1
    end
    
    methods
        function obj = untitled2(inputArg1,inputArg2)
            %UNTITLED2 Construct an instance of this class
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


classdef StackInterface
    %STACK interface abstract superclass to develop the concrete discrete
    %and continuous stacks.
    %   Subclasses should implement the properties and methods. The set/get
    %   access permission is still being decided.
    % Emilio Isaías-Camacho - Jan 2019 @GrohLabs
    
    properties
        Trigger
        ViewingWindow (1,2)
        Stack
    end
    
    methods
        function obj = StackInterface(inputArg1,inputArg2)
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


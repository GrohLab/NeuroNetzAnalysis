classdef GeneralWaveform
    %GENERALWAVEFORM implements a class for the
    %   Detailed explanation goes here
    
    properties 
        Data
        NSamples
        SamplingFreq
        Time
    end
    
    methods
        function obj = GeneralWaveform(data,samplingFreq)
            %UNTITLED11 Construct an instance of this class
            %   Detailed explanation goes here
                [rows, samples] = size(data);
                if rows > samples
                    data = data';
                end                
                obj.Data = data;
                obj.SamplingFreq = samplingFreq;
                
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


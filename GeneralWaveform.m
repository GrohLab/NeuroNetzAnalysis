classdef GeneralWaveform
    %GENERALWAVEFORM implements a class for the
    %   Detailed explanation goes here
    
    properties (SetAccess = 'private')
        Data 
        NSamples
        SamplingFreq
        Time
    end
    properties
        Units
    end
    
    methods
        function obj = GeneralWaveform(data, samplingFreq, units)
            %GENERALWAVEFORM Construct an instance of this class. It takes
            %two arguments: data and samplingfrequency which
            %   Detailed explanation goes here
                [rows, samples] = size(data);
                if rows > samples
                    data = data';
                end                
                obj.Data = data;
                obj.SamplingFreq = samplingFreq;
                obj.NSamples = length(data);
                obj.Time =...
                    0:1/obj.SamplingFreq:obj.(NSamples-1)/obj.SamplingFreq;
                if nargin == 3
                    if ischar(units)
                        obj.Units = units;
                    else
                        warning('String expected')
                    end
                end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        
    end
end


classdef StepWaveform < GeneralWaveform
    %STEPWAVEFORM is derived from the GeneralWaveform class and contains
    %only the extra necessary properties to produce time stamps for
    %triggering or anyother desired method.
    
    
    properties (Dependent)
        RiseAndFall
    end
    
    methods
        function obj = StepWaveform(data, samplingFreq, units, title)
            %STEPWAVEFORM Construct an instance of this class
            %   Detailed explanation goes here
            obj@GeneralWaveform(data,samplingFreq, units,title);
        end
        
        function  RaF = get.RiseAndFall(obj)
            ds = diff(obj.Data);
            rise = false(obj.NSamples,1);    % Rising edge times
            fall = rise;                    % Falling edge times
            % Maximum value divided by three
            rise(2:end) = ds > max(abs(ds))/3;
            fall(1:end-1) = ds < min(ds)/3;
            if sum(rise) ~= sum(fall)
                warning('The cardinality of the rising edges is different for the falling edges\n')
            else                
                RaF = [rise,fall];
            end
        end
    end
end


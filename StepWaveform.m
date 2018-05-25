classdef StepWaveform < DiscreteWaveform
    %STEPWAVEFORM is derived from the DiscreteWaveform class and contains
    %only the extra necessary properties to produce time stamps for
    %triggering or anyother desired method.
    
    properties
        TimeStamps
    end
    
    methods
        function obj = StepWaveform(data, samplingFreq, units, title)
            %STEPWAVEFORM Construct an instance of this class
            %   Detailed explanation goes here
            obj@DiscreteWaveform(data,samplingFreq, units,title);
        end
        function  RaF = get.TimeStamps(obj)
            ds = diff(obj.Data);
            rise = false(obj.NSamples,1);    % Rising edge times
            fall = rise;                    % Falling edge times
            % Maximum value divided by three
            rise(2:end) = ds > max(abs(ds))/3;
            rise = StepWaveform.CleanEdges(rise);
            fall(1:end-1) = ds < min(ds)/3;
            fall = StepWaveform.CleanEdges(fall);
            if sum(rise) ~= sum(fall)
                warning('The cardinality of the rising edges is different for the falling edges\n')
            else                
                RaF = [rise,fall];
            end
        end
    end
    methods (Static)
        function edgeOut = cleanEdges(edgeIn)
            doubleEdge = edgeIn(1:end-1) + edgeIn(2:end);
            repeatIdx = doubleEdge > 1;
            edgeOut = edgeIn;
            if sum(repeatIdx)
                edgeOut(repeatIdx) = false;
            end
        end
    end
end


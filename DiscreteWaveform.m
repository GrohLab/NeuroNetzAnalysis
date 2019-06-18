classdef (Abstract) DiscreteWaveform < GeneralWaveform
    %DISCRETEWAVEFORM This class would be abstract class to have
    %   Detailed explanation goes here
    
    properties (Abstract, SetAccess = 'private') 
        Triggers 
        FirstInTrain
        Count
        Delta
    end
    
    methods
        function obj = DiscreteWaveform(data, samplingFreq, units, title)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@GeneralWaveform(data, samplingFreq, units, title);
            
        end 
        function h = plot(obj,varargin)
            h = plot@GeneralWaveform(obj,varargin{:});
        end
    end
end

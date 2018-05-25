classdef (Abstract) DiscreteWaveform < GeneralWaveform
    %DISCRETEWAVEFORM This class would be abstract class to have
    %   Detailed explanation goes here
    
    properties (Abstract)
        TimeStamps
    end
    
    methods
        function obj = DiscreteWaveform(varargin)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@GeneralWaveform(varargin);
        end
        
    end
end

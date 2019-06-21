classdef (Abstract) DiscreteWaveform < GeneralWaveform
    %DISCRETEWAVEFORM This class would be abstract class to have
    %   Detailed explanation goes here
    
    properties (Abstract, SetAccess = 'private') 
        Triggers 
    end
    
    properties (SetAccess = 'private')
        FirstOfTrain
        LastOfTrain
        Count
        Delta
        
    end
    
    properties (Access = 'public')
        MinIEI (1,1) = [];
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
        
        % Get the first event of a train
        function fot = get.FirstOfTrain(obj)
            fs = obj.SamplingFreq;
            trigs = obj.Triggers;
            evntTms = trigs(:,1)./fs;
            if ~isempty(obj.MinIEI)
                fot = StepWaveform.firstOfTrain(evntTms,obj.MinIEI);
            else
                fprintf('The minimum inter-event interval MinIEI is not set!\n')
                fprintf('The function will use the mean MinIEI...\n')
                fot = StepWaveform.firstOfTrain(evntTms);
            end
        end
        
        % Get the last event of a train
        function lot = get.LastOfTrain(obj)
            fs = obj.SamplingFreq;
            trigs = obj.Triggers;
            evntTms = flip(trigs(:,1))./fs;
            if ~isempty(obj.MinIEI)
                lot = flip(StepWaveform.firstOfTrain(evntTms,obj.MinIEI));
            else
                fprintf('The minimum inter-event interval MinIEI is not set!\n')
                fprintf('The function will use the mean MinIEI...\n')
                lot = flip(StepWaveform.firstOfTrain(evntTms));
            end
        end
        
        % Get the number of events in a train
        function c = get.Count(obj)
            fot = obj.FirstOfTrain;
            lot = obj.LastOfTrain;
            feSubs = find(fot);
            leSubs = find(lot);
            c = (leSubs - feSubs) + 1;
        end
        
        % Get the duration of the train in seconds
        function d = get.Delta(obj)
            fot = obj.FirstOfTrain;
            lot = obj.LastOfTrain;
            feSubs = find(fot);
            leSubs = find(lot);
            d = abs(feSubs-leSubs)./obj.SamplingFreq;
        end
        
        
        
    end
end

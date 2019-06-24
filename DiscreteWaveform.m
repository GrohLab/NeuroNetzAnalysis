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
        MinIEI (1,1) = NaN; % Time in seconds
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
            if ~isnan(obj.MinIEI)
                fot = DiscreteWaveform.firstOfTrain(evntTms,obj.MinIEI);
            else
                fprintf('The minimum inter-event interval MinIEI is not set!\n')
                fprintf('The function will use the mean MinIEI...\n')
                [fot, obj.MinIEI] = DiscreteWaveform.firstOfTrain(evntTms);
            end
        end
        
        % Get the last event of a train
        function lot = get.LastOfTrain(obj)
            fs = obj.SamplingFreq;
            trigs = obj.Triggers;
            evntTms = flip(trigs(:,1))./fs;
            if ~isempty(obj.MinIEI)
                lot = flip(DiscreteWaveform.firstOfTrain(evntTms,obj.MinIEI));
            else
                fprintf('The minimum inter-event interval MinIEI is not set!\n')
                fprintf('The function will use the mean MinIEI...\n')
                [lot, obj.MinIEI] = flip(DiscreteWaveform.firstOfTrain(evntTms));
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
    
    %% Static public methods
    methods (Static, Access = 'public')
        % Recreate the signal with logic values.
        function logicalTrace = subs2idx(subs,N)
            if size(subs,2) == 2
                logicalTrace = false(1,N);
                for cmt = 1:size(subs,1)
                    logicalTrace(subs(cmt,1):subs(cmt,2)) = true;
                end
            else
                [Nr, Nc] = size(subs);
                logicalTrace = false(Nr * (Nr < Nc) + Nc * (Nc < Nr),...
                    N);
                subsClass = class(subs);
                switch subsClass
                    case 'cell'
                        for cmt = 1:size(subs,1)
                            logicalTrace(cmt,subs{cmt}) = true;
                        end
                    case {'single','double'}
                        logicalTrace(1,subs) = true;
                    otherwise
                        fprintf('Case not yet implemented...\n')
                end
            end
        end
        
        % Get a semi-logic trigger signal
        function semiLogicSignal = SCBrownMotion(RaF)
            [R,C] = size(RaF);
            if ~any([R == 2,C == 2])
                disp('What kind of rising and falling edges were given?')
                semiLogicSignal = RaF;
                return
            end
            if R > C
                RaF = RaF';
            end
            semiLogicSignal = cumsum(RaF(1,:) - RaF(2,:));
        end
        
        % Get the first true value of a logic pulse
        function [frstSpks, minIpi] = firstOfTrain(spkTimes, minIpi)
            % OUTPUTS a logical index for the edges which are the first in
            % time for the step pulse train.
            Ipi = abs(diff(spkTimes));
            if ~exist('minIpi','var')
                minIpi = mean(Ipi);
            end
            Pks = Ipi < minIpi;
            Sps = DiscreteWaveform.addFst(~Pks,true);
            frstSpks = DiscreteWaveform.addLst(Sps(1:end-1) & Pks,false);
        end 
    end
    
    %% Static private methods
    methods (Static, Access = 'private')
        function new_array = addFst(array,element)
            new_array = cat(find(size(array)~=1), element, array);
        end
        
        function new_array = addLst(array,element)
            new_array = cat(find(size(array)~=1), array, element);
        end
    end
end

classdef StepWaveform < DiscreteWaveform
    %STEPWAVEFORM is derived from the DiscreteWaveform class and contains
    %only the extra necessary properties to produce time stamps for
    %triggering or anyother desired method.
    
    properties
        Triggers
    end
    
    methods
        function obj = StepWaveform(data, samplingFreq, units, title)
            %STEPWAVEFORM Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                units = 'on-off';
                title = 'Step Waveform';
            end
            obj@DiscreteWaveform(data,samplingFreq, units,title);
            
        end
%         function set.Triggers(obj,trigs)
%             confirmation = input('Do you really want to overwrite the trigger times?','s');
%             if (confirmation == 'y' || confirmation == 'Y') &&...
%                     (isnumeric(trigs) && ~sum(sum(ceil(trigs) - floor(trigs))))
%                 obj.Triggers = trigs;
%             else
%                 disp('No changes made')
%             end
%         end
        function  RaF = get.Triggers(obj)
            if isa(obj.Data,'double')
                ds = diff(obj.Data);
                rise = false(obj.NSamples,1);    % Rising edge times
                fall = rise;                    % Falling edge times
                % Maximum value divided by three
                rise(2:end) = ds > max(abs(ds))/3;
                rise = StepWaveform.cleanEdges(rise);
                fall(1:end-1) = ds < min(ds)/3;
                fall = StepWaveform.cleanEdges(fall);
                if sum(rise) ~= sum(fall)
                    warning('The cardinality of the rising edges is different for the falling edges\n')
                else
                    RaF = [rise,fall];
                    obj.Triggers = RaF;
                end
            elseif isa(obj.Data,'logical')
                aux = obj.Data(1:end-1) - obj.Data(2:end);
                rise = find(aux == -1)' + 1;
                fall = find(aux == 1)';
                if (~isempty(rise) && isempty(fall)) || (numel(rise) ~= numel(fall))
                    fall = [fall;length(obj.Data)];
                end
                RaF = [rise,fall];
                obj.Triggers = RaF;
            end
        end
        function disp(obj)
            disp('Step waveform-------')
%             if ~isempty(obj.Triggers)
%                 fprintf('Title: %s\n',obj.Title)
%                 fprintf('Triggers: %d\n',length(obj.Triggers))
%                 fprintf('Sampling Frequency: %0.3f\n',obj.SamplingFreq/1e3)
%             end
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


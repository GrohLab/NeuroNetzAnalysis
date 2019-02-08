classdef StepWaveform < DiscreteWaveform
    %STEPWAVEFORM is derived from the DiscreteWaveform class and contains
    %only the extra necessary properties to produce time stamps for
    %triggering or anyother desired method.
    
    properties (SetAccess = 'private')
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
        
        function RaF = get.Triggers(obj)
            if isa(obj.Data,'double')
                ds = diff(obj.Data);
                if sum(ds>0) ~= numel(obj.Data)-1
                    rise = false(obj.NSamples,1);    % Rising edge times
                    fall = rise;                    % Falling edge times
                    % Maximum value divided by three
                    rise(2:end) = ds > max(abs(ds))/3;
                    rise = StepWaveform.cleanEdges(rise);
                    fall(1:end-1) = ds < min(ds)/3;
                    fall = StepWaveform.cleanEdges(fall);
                    if sum(rise) ~= sum(fall)
                        warning('The cardinality of the rising edges is different for the falling edges\n')
                        if numel(fall) < numel(rise)
                            fall = [fall;length(obj.Data)];
                        else
                            rise = [1;rise];
                        end
                    end
                    RaF = [rise,fall];
                    obj.Triggers = RaF;
                else
                    disp('The given data are already the triggers')
                    RaF = obj.Data;
                    obj.Triggers = RaF;
                end
            elseif isa(obj.Data,'logical')
                aux = obj.Data(1:end-1) - obj.Data(2:end);
                rise = find(aux == -1)' + 1;
                fall = find(aux == 1)';
                if (~isempty(rise) && isempty(fall)) || (numel(rise) ~= numel(fall))
                    if numel(fall) < numel(rise)
                        fall = [fall;length(obj.Data)];
                    else
                        rise = [1;rise];
                    end
                elseif sum((rise - fall) > 0)
                    if obj.Data(1)
                        rise = [1;rise];
                    end
                    if obj.Data(end)
                        fall = [fall;length(obj.Data)];
                    end
                elseif isempty(rise) && isempty(fall)
                    if sum(obj.Data)
                        % All data is a constant 1. Otherwise, the data
                        % contains no steps.
                        rise = 1;fall = obj.NSamples;
                    else
                        warning('No steps found in these data!')
                        fprintf('Returning empty variables...\n')
                    end
                end
                RaF = [rise,fall];
                obj.Triggers = RaF;
            end
        end
        function disp(obj)
            disp('Step waveform-------')
            if ~isempty(obj.Data)
                fprintf('Title: %s\n',obj.Title)
                fprintf('Triggers: %d\n',length(obj.Triggers))
                fprintf('Sampling Frequency: %0.3f kHz\n',obj.SamplingFreq/1e3)
            end
        end
        function myfunction(obj,inputArg1)
        end
    end
    methods (Static, Access = 'private')
        function edgeOut = cleanEdges(edgeIn)
            doubleEdge = edgeIn(1:end-1) + edgeIn(2:end);
            repeatIdx = doubleEdge > 1;
            edgeOut = edgeIn;
            if sum(repeatIdx)
                edgeOut(repeatIdx) = false;
            end
        end
    end
    methods (Static, Access = 'public')
        function logicalTrace = subs2idx(subs,N)
            if size(subs,2) == 2
                fprintf('Time windows\n')
                logicalTrace = false(1,N);
                for cmt = 1:size(subs,1)
                    logicalTrace(subs(cmt,1):subs(cmt,2)) = true;
                end
            else
                fprintf('Time points (i.e. spikes)\n')
                logicalTrace = false(size(subs,1),N);
                for cmt = 1:size(subs,1)
                    logicalTrace(cmt,subs{cmt}) = true;
                end
            end
        end
    end
end


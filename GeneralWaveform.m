classdef GeneralWaveform
    %GENERALWAVEFORM implements a class for the
    %   Detailed explanation goes here
    
    properties (SetAccess = 'private')
        Data 
        NSamples
        SamplingFreq
        Time
    end
    
    properties (Dependent, SetAccess = 'private')
        RiseAndFall
        Spikes
    end
    
    properties
        Units
        Title
    end
    
    
    methods
        function obj = GeneralWaveform(data, samplingFreq, units, title)
            %GENERALWAVEFORM Construct an instance of this class. It takes
            %two arguments: data and samplingfrequency which
            %   Detailed explanation goes here
            if nargin > 0 && ~isempty(data)
                [rows, samples] = size(data);
                if rows > samples
                    data = data';
                end                
                obj.Data = data;
                obj.SamplingFreq = samplingFreq;
                obj.NSamples = length(data);
                obj.Time =...
                    0:1/obj.SamplingFreq:(obj.NSamples-1)/obj.SamplingFreq;
                if nargin == 3
                    if ischar(units)
                        obj.Units = units;
                    else
                        warning('String expected; Units not assigned.')
                    end
                elseif nargin == 4
                    if ischar(title)
                        obj.Units = title;
                    else
                        warning('String expected; Title not assigned.')
                    end
                else
                    disp('Too many input arguments')
                end
            end
        end
        
        function obj = set.Title(obj,title)
            if ischar(title)
                obj.Title = title;
            else
                error('Expected char or string')
            end
        end
        
        function obj = set.Units(obj, units)
            if ischar(units)
                obj.Units = units;
            else
                error('Expected char or string')
            end
        end
        
        function h = plot(obj,varargin)
            %PLOT opens a new figure and plots the waveform with its
            %correct time and units.
            figure();h = plot(obj.Time,obj.Data,varargin{:});
            ylabel([obj.Units])
        end
        
        function RaF = get.RiseAndFall(obj)
            ds = diff(obj.Data);
            rise = false(obj.NSamples);    % Rising edge times
            fall = rise;                    % Falling edge times
            % Maximum value divided by three
            rise(2:end) = ds > max(abs(ds))/3;
            fall(1:end-1) = ds < min(ds)/3;
            if sum(rise) ~= sum(fall)
                warning('The cardinality of the rising edges is different for the falling edges\n')
            else
                RaF = cat(find(size(rise)==1),rise,fall);
                [Ns,~] = size(RaF);
                if Ns == 2
                    RaF = RaF';
                end
            end
        end
        function spkTimeStamps = getSpikesTimeStamps(obj,thresh,minISI)
            mx = 1/max(obj.Data);
            x=obj.Data*mx;
            dx=diff(x);
            if max(dx)>thresh
                %sp= x>thresh;
                %if ~sum(sp,'omitnan')
                sp = find(x>thresh);
                if ~isempty(sp)
                    if size(sp,1) < size(sp,2)
                        sp=sp';
                    end
                    dsp=[1000000;diff(sp)];
                    %indices=find(dsp>minISI);
                    indices=dsp>minISI;
                    
                    %indices
                    numel(sp);
                    sp=sp(indices);
                    starts=sp;
                    ends=sp+5;
                    
                    %ends(find(ends>numel(x)))=numel(x);
                    ends(ends>numel(x))=numel(x);
                    nsp=[];
                    for i=1:numel(sp)
                        v=x(starts(i):ends(i));
                        m=find(v==max(v));
                        m=m(1);
                        nsp(i)=starts(i)+m-1;
                    end
                    spkTimeStamps=nsp;
                end
            else
                spkTimeStamps=[];
            end
            obj.Spikes = spkTimeStamps;
        end
        
        function obj = set.Spikes(obj,spikes)
            obj.Spikes = spikes;
        end
    end
end


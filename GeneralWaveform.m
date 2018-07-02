classdef GeneralWaveform < handle
    %GENERALWAVEFORM implements a class for the any type of signals. This
    %class should convert to the discrete or continuous waveform the input
    %signals.
    %   Detailed explanation goes here
    
    properties (SetAccess = 'protected')
        SamplingFreq (1,1) double = 2e4;
        Data = [];
        NSamples
        Time
        ID
    end
    properties
        Units (1,:) char = 'mV';
        Title (1,:) char = '';
    end
    methods
        function obj = GeneralWaveform(data, samplingFreq, units, title)
            %GENERALWAVEFORM Construct an instance of this class. It takes
            %two mandatory arguments: data and samplingFreq, and two
            %optional: units and title, which can be modified afterwards.
            %   Detailed explanation goes here
            if nargin > 0 && ~isempty(data) && (isnumeric(data) || isa(data,'logical'))
                [rows, samples] = size(data);
                if rows > samples
                    data = data';
                end                
                if isnumeric(data)
                    obj.Data = double(data);
                elseif islogical(data)
                    obj.Data = logical(data);
                end
                obj.SamplingFreq = samplingFreq;
                obj.NSamples = length(data);
                obj.Time = seconds(...
                    0:1/obj.SamplingFreq:(obj.NSamples-1)/obj.SamplingFreq);
                if nargin >= 3
                    if ischar(units)
                        obj.Units = units;
                    else
                        warning('String expected; Units not assigned.')
                    end
                    if nargin == 4
                        if ischar(title)
                            obj.Title = title;
                        else
                            warning('String expected; Title not assigned.')
                        end
                    end
                elseif nargin > 4
                    disp('Too many input arguments')
                end
            else
                obj.Data = [];
            end
        end
        
        function set.Title(obj,title)
            if ischar(title)
                obj.Title = title;
            else
                error('Expected char or string')
            end
        end
        
        function set.Units(obj, units)
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
            ylabel(obj.Units);xlabel('Time (s)');title(obj.Title);
        end
        
        function set.SamplingFreq(obj,newFs)
            if newFs < 1
                warning('The given number is samller than 1. Assuming ''sampling period''')
                newFs = 1/newFs;
            end
            obj.SamplingFreq = newFs;
        end
        
    end
    methods (Abstract)
        disp(obj)
    end
end


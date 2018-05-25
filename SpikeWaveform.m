classdef SpikeWaveform < DiscreteWaveform
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess = 'private')
        % Multi or difficult spikes using Ultra Mega Sort 2000:
        UMSdata (1,1) UMSDataLoader;
        Bursts = BurstTrain();
    end
    properties
        TimeStamps
    end
    properties (Dependent)
        NumberOfSpikes;
    end
    properties (GetAccess = 'private',SetAccess = 'private')
        NormThreshold = [];
        MinISI (1,1) double = 1e-3;          % 1 ms
    end
    methods
        function obj = SpikeWaveform(varargin)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@DiscreteWaveform(varargin);
            if ~nargin
                disp('Using UMS interface')
            end
        end
        
        function obj = set.TimeStamps(obj, spks)
            obj.Spikes = spks;
        end
        
        % Function to use UMS to extract the spikes from the given data.
        function obj = getSpikes_UMS(varargin)
            obj.UMSdata = UMSDataLoader(varargin);
            obj.UMSdata.UMS2kPipeline;
            
        end
        function [spkTimeStamps, obj] = getSpikes_Thresh_ISI(obj,thresh,minISI)
            mx = 1/max(obj.Data);
            x=obj.Data*mx;
            disp(['Distance to the RMS: ',num2str(exp(-abs(rms(x)-thresh)))])
            % dx=diff(x);
            % if max(dx)>thresh % What does this if statement mean?
            if thresh < max(x)
                %sp= x>thresh;
                %if ~sum(sp,'omitnan')
                sp = find(x>thresh);
                if ~isempty(sp)
                    if size(sp,1) < size(sp,2)
                        sp=sp';
                    end
                    dsp=[inf;diff(sp)];
                    %indices=find(dsp>minISI);
                    indices=dsp>minISI;
                    
                    %indices
                    numel(sp);
                    sp=sp(indices);
                    starts=sp;
                    ends=sp+5;
                    
                    %ends(find(ends>numel(x)))=numel(x);
                    ends(ends>numel(x))=numel(x);
                    nsp=zeros(1,length(sp));
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
            obj.TimeStamps = spkTimeStamps;
        end
        
        function [Bursts, firstSpikes, C, BurstLength]=returnBursts(sp,isiCutoff)
            %[Bursts firstSpikes C]=returnBursts(sp,isiCutoff)
            %This function takes a spike train, sp, and an interspike interval
            %criterion, isiCutoff, and returns a structure of events, Bursts; each
            %entry represents a burst of spikes sorted according to isiCutoff.
            %firstSpikes gives a vector of the first spike in eache event.
            %C gives a vector of each event's size;
            %numel(C)=numel(firstSpikes)=numel(Bursts);
            BurstLength=[];
            if numel(sp)==0
                Bursts=[];
                firstSpikes=[];
                C=[];
            end
            % what if there is just one spike?
            if numel(sp)==1
                Bursts={sp}
                firstSpikes=sp;
                C=1;
            else
                %if there are multiple spikes, then sort them;
                dim=size(sp); if dim(2)>dim(1),sp=sp';end  %consistent column of spike times.
                Bursts={};
                %what if there are no spikes?
                if isempty(sp)
                    Bursts={};
                elseif numel(sp)==1   %this is redundant, but left for now.
                    Bursts{1}=sp;
                else
                    isis=[sp(1); diff(sp)];  %find interspike intervals
                    bursts=(find(isis>isiCutoff)); %indices of burst starts;
                    bursts=[1;bursts];bursts=unique(bursts);% what if first even is too close to beginning? count it as event.
                    firstSpikes=sp(bursts); %find first spikes;
                    C=nan(size(bursts)); %preallocate
                    for i=1:numel(bursts)
                        % for each burst starting spike, find spikes between it and the
                        % the next burst starting spike.
                        if i<numel(bursts)
                            i1=bursts(i);
                            i2=bursts(i+1)-1;
                            Bursts{i}=[sp(i1:i2)];
                        else
                            i1=bursts(i);
                            Bursts{i}=[sp(i1:end)];
                        end
                        C(i)=numel(Bursts{i});  %how many spikes per burst?
                    end
                end
            end
            if ~isempty(Bursts)
                for i=1:numel(Bursts)
                    BurstLength(i)=sum(diff(Bursts{i}));
                end
            end
            
        end
    end
end

function [expStack, LFPstack, Wstack] =...
    getStack(spT,alignP,ONOFF,timeSpan,fs,LFP,whiskerMovement,consEvents,fsLFP)
% GETSTACK returns a stack of spikes aligned to a certain event ''alignT''
% considering the events in the cell array ''consEvents''. The alignment
% can be done with the on-set or the off-set of the triggers using the
% ONOFF option. The default is 'on'. The timeSpan is given in a 1x2 array
% with the time before the trigger as the first element and the time after
% the trigger as the second element. The sampling frequency is always
% important and it is given in Hz. The LFP and whisker movement are
% examples of triggered averages from continuous signals. Finally, the
% consEvents are 

%% Computing the size of the PSTH stack
if isa(alignP,'logical')
    alWf = StepWaveform(alignP,fs,'on-off','Align triggers');
    alignP = alWf.Triggers;
end
switch ONOFF
    case 'on'
        disp('Considering onset of the triggers')
    case 'off'
        disp('Considering offset of the triggers')
    otherwise
        disp('Unrecognized trigger selection. Considering onsets')
        ONOFF = 'on';
end
[auxR,auxC] = size(alignP);
if auxR < auxC
    alignP = alignP';
end
[Na, raf] = size(alignP);
if raf > 2 || raf < 1
    fprintf(['Warning! The alignment matrix is expected to have ',...
        'either only rising or rising and falling edges time indices.\n'])
    expStack = NaN;
    return;
end

%% Considered events arrangement.
Ne = 0;
if exist('consEvents','var') && ~isempty(consEvents)
    typ = whos('consEvents');
    switch typ.class
        case 'double'
            [~, Ne] = size(consEvents);
            if mod(Ne,2)
                Ne = Ne/2;
            else
                fprintf('Omitting the events to consider.\n')
                Ne = 0;
            end
        case 'cell'
            Ne = length(consEvents);
            consEvents2 = consEvents;
            evntTrain = cellfun(@islogical,consEvents);
            % Converting the logical event trains into indices
            for ce = 1:Ne
                if evntTrain(ce)
                    stWv = StepWaveform(consEvents{ce},fs);
                    consEvents2{ce} = stWv.Triggers;
                end
            end
            consEvents = consEvents2;
        otherwise
            fprintf('The events to consider are not in recognized format.\n')
    end
end
%% Preallocation of the spike-stack:
toi = sum(timeSpan);
prevSamples = ceil(timeSpan(1) * fs);
postSamples = ceil(timeSpan(2) * fs);
%%%%%%%%%%%%%%%%%%%%%%%%%%% BEWARE OF THE LFP SAMPLING FREQUENCY currently
%%%%%%%%%%%%%%%%%%%%%%%%%%% set at 1000 Hz
% fsLFP = 1e3;
fsConv = fsLFP/fs;
%Nt = round(toi*fs)+1;
Nt = prevSamples + postSamples + 1;
expStack = false(2+Ne,Nt,Na);
prevSamplesLFP = ceil(timeSpan(1) * fsLFP);
postSamplesLFP = ceil(timeSpan(2) * fsLFP);
NtLFP = round(toi*fsLFP) + 1;
LFPstack = zeros(NtLFP,Na);
Wstack = LFPstack;
if isnumeric(spT)
    mxS = spT(end) + Nt;
    spTemp = false(1,mxS);
    spTemp(spT) = true;
    spT = spTemp;
end
%% Cutting the events into the desired segments.
for cap = 1:Na
    % Considering the rising or the falling edge of the step function.
    if strcmp(ONOFF, 'on')
        segmIdxs = [alignP(cap,1)-prevSamples,alignP(cap,1)+postSamples];
        segmIdxsLFP = round([(alignP(cap,1)*fsConv)-prevSamplesLFP,...
            (alignP(cap,1)*fsConv)+postSamplesLFP]);
    elseif strcmp(ONOFF, 'off')
        segmIdxs = [alignP(cap,2)-prevSamples,alignP(cap,2)+postSamples];
        segmIdxsLFP = round([(alignP(cap,2)*fsConv)-prevSamplesLFP,...
            (alignP(cap,2)*fsConv)+postSamplesLFP]);
    end
    % The segments should be in the range of the spike train.
    if segmIdxs(1) >= 1 && segmIdxs(2) <= length(spT)
        spSeg = spT(segmIdxs(1):segmIdxs(2));
        if ~(segmIdxsLFP(1) >= 1 && segmIdxsLFP(2) <= length(LFP)) ||...
                ~(segmIdxsLFP(1) >= 1 && segmIdxsLFP(2) <= length(whiskerMovement))
            OUTFLAG = true;
        else
            OUTFLAG = false;
        end
    else
        Na = Na - 1;
        continue;
    end
    expStack(2,:,cap) = spSeg;
    % Find 'overlapping' periods in time of interest
    alignPeriod = getEventPeriod(alignP, {alignP}, ONOFF, cap,...
        prevSamples, postSamples);
    expStack(1,:,cap) = alignPeriod;
    if Ne
        expStack(3:2+Ne,:,cap) =...
            getEventPeriod(alignP, consEvents, ONOFF, cap,...
            prevSamples, postSamples);
    end
    % Getting the LFP segments taking into account the different sampling
    % frequencies i.e. LFP-->1 kHz Spikes --> 20 kHz HARD CODE!! BEWARE!!
    if ~OUTFLAG
        if exist('LFP','var') && ~isempty(LFP)
            LFPstack(:,cap) = LFP(round((segmIdxsLFP(1):segmIdxsLFP(2))));
        end
        if exist('whiskerMovement','var') && ~isempty(whiskerMovement)
            Wstack(:,cap) = whiskerMovement(round((segmIdxsLFP(1):segmIdxsLFP(2))));
        end
    else
        disp('What to do in case that the indexes are outside the window?')
    end
end
end

% Aligning the events according to the considered time point. The inputs
% are time indices called Tdx as in Time inDeX.
function evntOn = getEventPeriod(alignTdx, evntTdx, ONOFF, cap, prev, post)
% Assuming that the considered events are always cells.
if isempty(evntTdx)
    evntOn = [];
    return;
else
    evntOn = false(numel(evntTdx),prev+post+1);
    for ce = 1:length(evntTdx)
        % Change of implementation. Considering the on set OR the offset of
        % a step pulse.
        if strcmp(ONOFF,'on')
            relTdx = evntTdx{ce} - alignTdx(cap,1);
        else
            relTdx = evntTdx{ce} - alignTdx(cap,2);
        end
        
        % Indexes for encountered rising (onIdx) and falling (offIdx) edges
        onIdx = relTdx(:,1) >= -prev & relTdx(:,1) < post; 
        if size(relTdx, 2) == 2
            offIdx = relTdx(:,2) >= -prev & relTdx(:,2) < post;
        else
            offIdx = onIdx;
        end
        % Indexes considering the partial or complete step that falls into
        % the considered segment.
        allIdx = find(onIdx | offIdx);
        initStep = relTdx(onIdx,1) + prev + 1;
        fnalStep = relTdx(offIdx,2) + prev + 1;
        % Event rising, falling or both into the alignment window as a
        % first condition. 
        stCount = 1;
        ndCount = 1;
        for cstp = 1:numel(allIdx)
            % If there is no rising edge found, then the step is considered
            % to start before the window and it will be true since the
            % start of the window. Otherwise, the obvious index is taken.
            if onIdx(allIdx(cstp))
                strt = initStep(stCount);
                stCount = stCount + 1;
            else
                strt = 1;
            end
            % Same case. If the falling edge is not found, the step is
            % considered to end in a later time point than the considered
            % window and it will be true from the rising edge to the end.
            if offIdx(allIdx(cstp))
                fnal = fnalStep(ndCount);
                ndCount = ndCount + 1;
            else
                fnal = post + prev + 1;
            end
            evntOn(ce,strt:fnal) = true;
        end
    end
end
end


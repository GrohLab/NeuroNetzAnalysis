function [discreteStack, continuouStack] =...
    getStacks(spT,alignP,ONOFF,timeSpan,fs,fsLFP,consEvents,varargin)
% GETSTACK returns a stack of spikes aligned to a certain event ''alignT''
% considering the events in the cell array ''consEvents''. The alignment
% can be done with the on-set or the off-set of the triggers using the
% ONOFF option. The default is 'on'. The timeSpan is given in a 1x2 array
% with the time before the trigger as the first element and the time after
% the trigger as the second element. The sampling frequency is always
% important and it is given in Hz. The LFP and whisker movement are
% examples of triggered averages from continuous signals. Finally, the
% consEvents are 
discreteStack = NaN;
continuouStack = NaN;
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
% if auxR < auxC
%     alignP = alignP';
% end
Na = auxR;
if auxC > 2 || auxC < 1
    fprintf(['Warning! The alignment matrix is expected to have ',...
        'either only rising or rising and falling edges time indices.\n'])
    return;
end

%% Considered events arrangement.
Ne = 0;
if exist('consEvents','var') && ~isempty(consEvents)
    typ = whos('consEvents');
    switch typ.class
        case {'double', 'single'}
            [rws, cols] = size(consEvents);
            if rws < cols
                Ne = rws;
                consEvents = consEvents';
            else
                Ne = cols;
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
%% Preallocation of the discrete stack:
prevSamples = ceil(timeSpan(1) * fs);
postSamples = ceil(timeSpan(2) * fs);
Nt = prevSamples + postSamples + 1;
discreteStack = false(2+Ne,Nt,Na);
% Creation of the logical spike train
if isnumeric(spT)
    mxS = spT(end) + Nt;
    spTemp = false(1,mxS);
    spTemp(spT) = true;
    spT = spTemp;
end
%% Preallocation of the continuous stack:
if ~exist('fsLFP','var')
    fsLFP = fs;
end 
fsConv = fsLFP/fs;
% Signal validation
Ns = numel(varargin);
if Ns
    signalCheck = cellfun(@isnumeric,varargin);
    signalCheck2 = cellfun(@length,varargin);
    if sum(signalCheck) ~= Ns
        fprintf('Discarding those inputs which are not numeric...\n')
        disp(varargin(~signalCheck))
        varargin(~signalCheck) = [];
        Ns = Ns - sum(~signalCheck);
        signalCheck2(~signalCheck) = [];
    end
    if std(signalCheck2) ~= 0
        fprintf('The signals have not the same length...\n')
        fprintf('Considering the smallest: %d\n',min(signalCheck2))
        MAX_CONT_SAMP = min(signalCheck2);
    else
        MAX_CONT_SAMP = signalCheck2(1);
    end
    prevSamplesLFP = ceil(timeSpan(1) * fsLFP);
    postSamplesLFP = ceil(timeSpan(2) * fsLFP);
    NtLFP = prevSamplesLFP + postSamplesLFP + 1;
    continuouStack = single(zeros(Ns,NtLFP,Na));
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
    % Validations for the discrete stack
    spSeg = false(1,Nt);
    if segmIdxs(1) >= 1 && segmIdxs(2) <= length(spT)
        spSeg = spT(segmIdxs(1):segmIdxs(2));
    elseif segmIdxs(1) <= 0
        fprintf('The viewing window is out of the signal range.\n')
        fprintf('%d samples required before the signal. Omitting...\n',...
            segmIdxs(1))
        segmIdxs(segmIdxs <= 0) = 1;
        spSeg(Nt - segmIdxs(2) + 1:Nt) = spT(segmIdxs(1):segmIdxs(2));
    else
        fprintf('The viewing window is out of the signal range.\n')
        fprintf('%d samples required after the signal. Omitting...\n',...
            segmIdxs(2))
        segmIdxs(segmIdxs > length(spT)) = length(spT);
        spSeg(1:Nt - diff(segmIdxs)) = spT(segmIdxs(1):segmIdxs(2));
    end
    discreteStack(2,:,cap) = spSeg;
    % Find 'overlapping' events in time of interest
    alignPeriod = getEventPeriod(alignP, {alignP}, ONOFF, cap,...
        prevSamples, postSamples);
    discreteStack(1,:,cap) = alignPeriod;
    % If there are continuous events given
    if Ne
        if isa(consEvents,'cell')
        discreteStack(3:2+Ne,:,cap) =...
            getEventPeriod(alignP, consEvents, ONOFF, cap,...
            prevSamples, postSamples);
        elseif isnumeric(consEvents)
            discreteStack(3:2+Ne,:,cap) =...
            getEventPeriod(alignP, {consEvents}, ONOFF, cap,...
            prevSamples, postSamples);
        else
        end
    end
    % Getting the continuous segments if the time window is in range of the
    % signal domain and validation for the continuous stack
    if Ns
        % segmIdxsLFP(1) >= 1 && segmIdxsLFP(2) <= MAX_CONT_SAMP
        
        % continuouStack(:,:,cap) = single(cell2mat(signalSegments'));
        continuouStack(:,:,cap) =...
            getContinuousSignalSegment(varargin, segmIdxsLFP, MAX_CONT_SAMP, NtLFP);
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
        
        % Indexes for encountered rising _onIdx_, falling _offIdx_ edges,
        % or a bigger pulse envolving the viewing window _invIdx_.
        onIdx = relTdx(:,1) >= -prev & relTdx(:,1) < post; 
        if size(relTdx, 2) == 2
            offIdx = relTdx(:,2) >= -prev & relTdx(:,2) < post;
            invIdx = relTdx(:,1) < -prev & relTdx(:,2) > post;
        else
            offIdx = onIdx;
            invIdx = onIdx;
        end
        % Indexes considering the partial or complete step that falls into
        % the considered segment.
        allIdx = find(onIdx | offIdx | invIdx);
        initStep = relTdx(onIdx,1) + prev + 1;
        if size(relTdx, 2) == 2
            fnalStep = relTdx(offIdx,2) + prev + 1;
        else
            fnalStep = initStep;
        end
        % Indices considering the inclusion of the considered window into a
        % lengthy pulse of the considered event.
        %
        % HERE
        %
        
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

function contSigSeg = getContinuousSignalSegment(signalCell, Idxs, N, Nt)
contSigSeg = zeros(numel(signalCell),Nt,'single');
if Idxs(1) >= 1 && Idxs(2) <= N
    signalSegments = getSignalSegments(signalCell, Idxs);
    contSigSeg = single(cell2mat(signalSegments'));
elseif Idxs(1) <= 0
    Idxs(Idxs <= 0) = 1;
    signalSegments = getSignalSegments(signalCell, Idxs);
    contSigSeg(:,Nt - Idxs(2) + 1:Nt) = single(cell2mat(signalSegments'));
else
    Idxs(Idxs > N) = N;
    signalSegments = getSignalSegments(signalCell, Idxs);
    contSigSeg(:,1:diff(Idxs)+1) = single(cell2mat(signalSegments'));
end
end

function sigSeg = getSignalSegments(signalCell,Idxs)
try
    sigSeg =...
        cellfun(...
        @(x) x(round((Idxs(1):Idxs(2)))),...
        signalCell, 'UniformOutput', false);
catch ME
    disp(ME.getReport)
    fprintf('Very unlikely case: signals with different length\n')
    fprintf('Worth debugging!\n')
end
transpSign = cellfun(@isrow,sigSeg);
if sum(~transpSign)
    try
        sigSeg(~transpSign) =...
            {sigSeg{~transpSign}'};
    catch
        auxSegm = cellfun(@transpose,sigSeg(~transpSign),...
            'UniformOutput',false);
        sigSeg(~transpSign) = auxSegm;
    end
end
end


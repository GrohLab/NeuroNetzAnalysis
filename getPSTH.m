function [PSTH, figID] = getPSTH(PSTHstack, plotQ, saveQ, cellType, binSz, fs)
% GETPSTH returns a peri-stimulus triggered histogram given a triggered
% aligned stack in the form MxNxT, where M is the number of channels to
% align, N is the number of samples (time) and T is the number of aligning
% triggers found in the channel. 
% The kicking out process needs to be implemented in a dynamical form. The
% events after the (first) spike channel are a guide to kick row of the
% spikes out. For now, we will kick out all of the rows which contain a
% true in the observed time. 

% Computing the size of the PSTH stack
if nargin < 6
    fs = 2e4;
    if nargin < 5
        binSz = 1/2e4;
        if nargin < 4
            cellType = 'POm';
            if nargin < 3
                saveQ = false;
                if nargin < 2
                    plotQ = false;
                end
            end
        end
    end
end
[Ne, Nt, Na] = size(PSTHstack);

% Kicking out everything that contains a true in the channel. No light, no
% puff, no touch. Nothing.
kickOutIdx = [];
if size(PSTHstack,1) > 2
    kickOutIdx = squeeze(sum(PSTHstack(3:end,:,:),2));
    kickOutIdx(kickOutIdx ~= 0) = true;
end
PSTH = zeros(2,ceil(Nt/binSz));
PSTH = squeeze(sum(PSTHstack(:,kickOutIdx),3));
binEls = round(binSz * fs);
cb = 0;
while cb < Naux/binEls - 1
    PSTH(cb+1) = sum(auxCounts(cb*binEls+1:(cb+1)*binEls));
    cb = cb + 1;
end
if plotQ
    switch cellType
        case 'POm'
            % Black for POm
            colr = 'k';
        case 'VPM'
            % Gray for VPM
            colr = 0.6 * ones(1,3);
        otherwise
            % Cyan for other.
            colr = [77, 210, 255]/255;
    end
    % Create figure for the PSTH. The color is determined by the cell type
    figID = figure('Name',[cellType, ' PSTH'],'Color',[1,1,1]);
    
else
    return;
end


[Na, raf] = size(alignP);
if raf > 2 || raf < 1
    fprintf(['Warning! The alignment matrix is expected to have ',...
        'either only rise or rise and falling edges time indices.\n'])
    PSTHstack = NaN;
    return;
end
Ne = 0;
if nargin == 6
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
% Preallocation of the spike-stack:
toi = sum(timeSpan);
prevSamples = ceil(timeSpan(1) * fs);
postSamples = ceil(timeSpan(2) * fs);
Nt = round(toi*fs) + 1;
PSTHstack = false(2+Ne,Nt,Na);
LFPstack = zeros(Nt,Na);
Wstack = LFPstack;
if isnumeric(spT)
    mxS = spT(end) + Nt;
    spTemp = false(1,mxS);
    spTemp(spT) = true;
    spT = spTemp;
end
for cap = 1:Na
    segmIdxs = [alignP(cap,1)-prevSamples,alignP(cap,1)+postSamples];
    % The segments should be in the range of the spike train.
    if segmIdxs(1) >= 1 && segmIdxs(2) <= length(spT)
        spSeg = spT(segmIdxs(1):segmIdxs(2));
    else
        Na = Na - 1;
        continue;
    end
    PSTHstack(2,:,cap) = spSeg;
    % Find 'overlapping' periods in time of interest
    alignPeriod = getEventPeriod(alignP,{alignP},cap,prevSamples,postSamples);
    PSTHstack(1,:,cap) = alignPeriod;
    if Ne
        PSTHstack(3:2+Ne,:,cap) =...
            getEventPeriod(alignP, consEvents, cap, prevSamples, postSamples);
    end
    % Getting the LFP segments taking into account the different sampling
    % frequencies i.e. LFP-->1 kHz Spikes --> 20 kHz HARD CODE!! BEWARE!!
    if exist('LFP','var') && ~isempty(LFP)
        LFPstack(:,cap) = LFP(round((segmIdxs(1):segmIdxs(2))*(1e3/fs)));
    end
    if exist('whiskerMovement','var') && ~isempty(whiskerMovement)
        Wstack(:,cap) = whiskerMovement(round((segmIdxs(1):segmIdxs(2))*(1e3/fs)));
    end
end
end

% Aligning the events according to the considered time point. The inputs
% are time indices called Tdx as in Time inDeX.
function evntOn = getEventPeriod(alignTdx, evntTdx, cap, prev, post)
% Assuming that the considered events are always cells.
if isempty(evntTdx)
    evntOn = [];
    return;
else
    evntOn = false(numel(evntTdx),prev+post+1);
    for ce = 1:length(evntTdx)
        relTdx = evntTdx{ce}(:,1) - alignTdx(cap,1);
        inToi = find(relTdx >= -prev & relTdx < post);
        if ~isempty(inToi)
            psthIdx = evntTdx{ce}(inToi,1)-alignTdx(cap,1) + prev + 1;
            lenIdx = evntTdx{ce}(inToi,2);
            for cep = 1:numel(psthIdx)
                idxs = psthIdx(cep):psthIdx(cep)+lenIdx(cep)-1;
                idxs = idxs(idxs <= post+prev+1);
                evntOn(ce,idxs) = true;
            end
        end
    end
end
end


% function [PSTH,Trials] = getPSTH(spT,stimPeriods,whp,puffPeriods,LFP,...
%     winPre,winPost,fs)
% % GETSTA gets an average of the stimulus $s$ starting at time $t$ and
% % ending at time $t_0$, with a sampling period of $\delta$T. The time point
% % $t_0$ corresponds to the considered spike time $sp_i$ and the absolute
% % difference between $t$ and $t_0$ is 'win' in \milli seconds. Usually 10
% % ms.
% winel = round(winPre*fs); % window elements
% postSamples = round(winPost*fs);
% len = length(LFP);
% [PSTH,Trials] = getPeriStimuliSpikes(...
%     winel,postSamples,stimPeriods,whp,puffPeriods,spT,len,fs);
%
% end
%
% function [SpSum, Trials] = getPeriStimuliSpikes(prevSamples,...
%     postSamples,stimTime,wh,puffPeriods,logSpTrain,sigLen,fs)
% % Initializing the output variables.
% puffOmit = 0.2; % Puff omitting window 100 ms
% SpSumW = zeros(1,prevSamples+postSamples+1);
% SpSumNW = zeros(1,prevSamples+postSamples+1);
% SpSumA = zeros(1,prevSamples+postSamples+1);
% Nst = sum(stimTime(1,:)>0);
% WT = 0;
% NWT = 0;
% AT = 0;
% puffLight = 0;
% for cst = 1:Nst
%     % Taking a look only at light stimuli lasting 0.5 seconds
%     if stimTime(2,cst) > 0.4 %&& stimTime(2,cst) < 0.6
%         segmIdxs = [stimTime(1,cst)-prevSamples,stimTime(1,cst)+postSamples];
%         if segmIdxs(1) < 1 || segmIdxs(2) > sigLen
%             fprintf('Window is outside the signal domain. Ignoring...\n')
%             continue;
%         else
%             if ~sum(puffPeriods(stimTime(1,cst)-round(puffOmit*fs):...
%                     stimTime(1,cst)+round(puffOmit*fs)))
%                 spTrainSegment = logSpTrain(segmIdxs(1):segmIdxs(2));
%                 SpSumA = SpSumA + spTrainSegment;
%                 AT = AT + 1;
%                 if wh(stimTime(1,cst))
%                     WT = WT + 1;
%                     SpSumW = SpSumW + spTrainSegment;
%                 else
%                     NWT = NWT + 1;
%                     SpSumNW = SpSumNW + spTrainSegment;
%                 end
%             else
%                 puffLight = puffLight + 1;
%             end
%         end
%     end
% end
% fprintf('Puff: %d\n',puffLight)
% SpSum = [SpSumA;SpSumW;SpSumNW];
% Trials = [AT;WT;NWT];
% end

% function [pB,pT,pS] = estimateStimuliPDF(burstStim,tonicStim,stim)
% M = 3;
% err = 1e-3;
% xLow = min(stim);
% xHigh = max(stim);
% dataRange = [xLow,xHigh];
% prB = emforgmm(burstStim,M,err,0);        % Stimuli PDF for busrts
% prT = emforgmm(tonicStim,M,err,0);  % Stimuli PDF for tonic sp
% prS = emforgmm(stim,M,err,0);        % Stimuli PDF for non-sp
% pB = gmmpdf(prB,dataRange);
% pT = gmmpdf(prT,dataRange);
% pS = gmmpdf(prS,dataRange);
% end
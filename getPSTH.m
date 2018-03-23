function [PSTHstack, LFPstack, Wstack] =...
    getPSTH(spT,alignP,timeSpan,fs,consEvents)
% GETPSTH returns a stack of spikes aligned to a certain event ''alignT''
% considering the events in the cell array ''consEvents''.

% Computing the size of the PSTH stack
[auxR,auxC] = size(alignP);
if auxR < auxC
    alignP = alignP';
end
[Na, raf] = size(alignP);
if raf > 2 || raf < 1
    fprintf(['Warning! The alignment matrix is expected to have ',...
        'either only rise or rise and falling edges time indices.\n'])
    PSTHstack = NaN;
    return;
end
Ne = 0;
if nargin == 5
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
function [PSTH,WT,NWT] = getPSTH(spT,stimPeriods,whp,LFP,winPre,winPost,fs)
% GETSTA gets an average of the stimulus $s$ starting at time $t$ and
% ending at time $t_0$, with a sampling period of $\delta$T. The time point
% $t_0$ corresponds to the considered spike time $sp_i$ and the absolute
% difference between $t$ and $t_0$ is 'win' in \milli seconds. Usually 10
% ms.
winel = round(winPre*fs); % window elements
postSamples = round(winPost*fs);
len = length(LFP);
[PSTH,WT,NWT] = getPeriStimuliSpikes(winel,postSamples,stimPeriods,whp,spT,fs,len);

end

function [SpSum, WT, NWT] = getPeriStimuliSpikes(prevSamples,...
    postSamples,stimTime,wh,spT,fs,sigLen)
SpSumW = zeros(1,prevSamples+postSamples+1);
SpSumNW = zeros(1,prevSamples+postSamples+1);
Nst = sum(stimTime>0);
logSpTrain = false(1,sigLen);
logSpTrain(round(spT*fs)) = true;
WT = 0;
NWT = 0;
for cst = 1:Nst
    segmIdxs = [stimTime(cst)-prevSamples,stimTime(cst)+postSamples];
    if segmIdxs(1) < 1 || segmIdxs(2) > sigLen
        fprintf('hehe the window is outside the signal domain. We will just ignore it.\n')
        continue;
    else
        spTrainSegment = logSpTrain(segmIdxs(1):segmIdxs(2));
        SpSumW = SpSumW + spTrainSegment;
        WT = WT + 1;
%         if wh(stimTime(cst))
%             WT = WT + 1;
%             SpSumW = SpSumW + spTrainSegment;
%         else
%             NWT = NWT + 1;
%             SpSumNW = SpSumNW + spTrainSegment;
%         end
    end
end
SpSum = [SpSumW;SpSumNW];
end

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
function [PSTH,Trials] = getPSTH(spT,stimPeriods,whp,puffPeriods,LFP,...
    winPre,winPost,fs)
% GETSTA gets an average of the stimulus $s$ starting at time $t$ and
% ending at time $t_0$, with a sampling period of $\delta$T. The time point
% $t_0$ corresponds to the considered spike time $sp_i$ and the absolute
% difference between $t$ and $t_0$ is 'win' in \milli seconds. Usually 10
% ms.
winel = round(winPre*fs); % window elements
postSamples = round(winPost*fs);
len = length(LFP);
[PSTH,Trials] = getPeriStimuliSpikes(...
    winel,postSamples,stimPeriods,whp,puffPeriods,spT,len,fs);

end

function [SpSum, Trials] = getPeriStimuliSpikes(prevSamples,...
    postSamples,stimTime,wh,puffPeriods,logSpTrain,sigLen,fs)
% Initializing the output variables.
puffOmit = 0.2; % Puff omitting window 100 ms
SpSumW = zeros(1,prevSamples+postSamples+1);
SpSumNW = zeros(1,prevSamples+postSamples+1);
SpSumA = zeros(1,prevSamples+postSamples+1);
Nst = sum(stimTime(1,:)>0);
WT = 0;
NWT = 0;
AT = 0;
puffLight = 0;
for cst = 1:Nst
    % Taking a look only at light stimuli lasting 0.5 seconds
    if stimTime(2,cst) > 0.4 %&& stimTime(2,cst) < 0.6
        segmIdxs = [stimTime(1,cst)-prevSamples,stimTime(1,cst)+postSamples];
        if segmIdxs(1) < 1 || segmIdxs(2) > sigLen
            fprintf('Window is outside the signal domain. Ignoring...\n')
            continue;
        else
            if ~sum(puffPeriods(stimTime(1,cst)-round(puffOmit*fs):...
                    stimTime(1,cst)+round(puffOmit*fs)))
                spTrainSegment = logSpTrain(segmIdxs(1):segmIdxs(2));
                SpSumA = SpSumA + spTrainSegment;
                AT = AT + 1;
                if wh(stimTime(1,cst))
                    WT = WT + 1;
                    SpSumW = SpSumW + spTrainSegment;
                else
                    NWT = NWT + 1;
                    SpSumNW = SpSumNW + spTrainSegment;
                end
            else
                puffLight = puffLight + 1;
            end
        end
    end
end
fprintf('Puff: %d\n',puffLight)
SpSum = [SpSumA;SpSumW;SpSumNW];
Trials = [AT;WT;NWT];
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
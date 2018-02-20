function [STA, STSTD] = getSTA(sp,s,maxISI,win,fs)
% GETSTA gets an average of the stimulus $s$ starting at time $t$ and
% ending at time $t_0$, with a sampling period of $\delta$T. The time point
% $t_0$ corresponds to the considered spike time $sp_i$ and the absolute
% difference between $t$ and $t_0$ is 'win' in \milli seconds. Usually 10
% ms.
[row,col] = size(s);
if row > col
    s = s';
end
[bT, tT] = getInitialBurstSpike(sp,maxISI*fs);
spT = [bT,tT];          % Stack first bursts and then tonic spikes
Ns = length(spT);
Nb = length(bT);
winel = round(win*fs); % window elements
stimStack = zeros(Ns,winel);
stimNotSpike = zeros(1,Ns*winel);
suma = 1;
for csp = 1:Ns
    stimStack(csp,:) = s(spT(csp)-winel+1:spT(csp));
    if csp == 1
        aux = s(1:spT(csp)-winel);
    else
        aux = s(spT(csp-1)+1:spT(csp)-winel);
    end
    lngth = length(aux);
    stimNotSpike(suma:suma+lngth-1) = aux;
    suma = suma + lngth;
end
[pB, pT, pS]=...       Bursts            Tonic Spikes          Non-spiking
    estimateStimuliPDF(stimStack(1:Nb,:),stimStack(Nb+1:end,:),stimNotSpike);
KL_SB = KullbackLeiblerDivergence(pS,pB);
KL_ST = KullbackLeiblerDivergence(pS,pT);
STAb = mean(stimStack(1:Nb,:),1);               % Burst spikes STA
STAt = mean(stimStack(Nb+1:end,:),1);           % Tonic spikes STA
STSTDb = std(stimStack(1:Nb,:),[],1);
STSTDt = std(stimStack(Nb+1:end,:),[],1);
end

function [pB,pT,pS] = estimateStimuliPDF(burstStim,tonicStim,stim)
M = 3;
err = 1e-3;
prB = emforgmm(burstStim,M,err,0);        % Stimuli PDF for busrts
prT = emforgmm(tonicStim,M,err,0);  % Stimuli PDF for tonic sp
prS = emforgmm(stim,M,err,0);        % Stimuli PDF for non-sp
pB = getPDFFromGMM(prB,burstStim);
pT = getPDFFromGMM(prT,tonicStim);
pS = getPDFFromGMM(prS,stim);
end

function plotSTA(STA, stimStack, STSTD)

end

function [bIdx, tIdx] = getInitialBurstSpike(spT,maxISI)
% GETINITIALBUSRTSPIKE gets the first spike time for each burst and the
% tonic spikes are unmodified according to the maximum
% inter-spiking-interval.
delta_spT = [inf,diff(spT)];
fsIdx = delta_spT >= maxISI;    % First spike index (burst-wise)
bsIdx = [~fsIdx(2:end) & fsIdx(1:end-1),fsIdx(end)];
tsIdx = ~bsIdx & fsIdx;
bIdx = spT(bsIdx);
tIdx = spT(tsIdx);
end
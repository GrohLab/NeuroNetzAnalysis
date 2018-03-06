function [STA, STSTD] = getSTA(sp,s,maxISI,win,fs)
% GETSTA gets an average of the stimulus $s$ starting at time $t$ and
% ending at time $t_0$, with a sampling period of $\delta$T. The time point
% $t_0$ corresponds to the considered spike time $sp_i$ and the absolute
% difference between $t$ and $t_0$ is 'win' in \milli seconds. Usually 10
% ms.

% NEED TO TAKE ALL THE INPUTS INTO CONSIDERATION!!

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
stimIdxs = false(size(s));
for csp = 1:Ns
    stimStack(csp,:) = s(spT(csp)-winel+1:spT(csp));
    stimIdxs(spT(csp)-winel+1:spT(csp)) = true;
    disp(['Window size: ',sprintf('%d - %d',spT(csp)-winel+1,spT(csp))])
    if csp == 1
        aux = s(1:spT(csp)-winel);
    else
        aux = s(spT(csp-1)+1:spT(csp)-winel);
    end
    lngth = length(aux);
    stimNotSpike(suma:suma+lngth-1) = aux;
    suma = suma + lngth;
end
%% Gaussian Mixture Model estimation
% [pB, pT, pS]=...       Bursts            Tonic Spikes          Non-spiking
%     estimateStimuliPDF(stimStack(1:Nb,:),stimStack(Nb+1:end,:),stimNotSpike);
% KL_SB = pS.comparePDF(pB);
% KL_ST = pS.comparePDF(pT);
% Use the histogram values instead of estimating a pdf.
%%
STAb = mean(stimStack(1:Nb,:),1);               % Burst spikes STA
STAt = mean(stimStack(Nb+1:end,:),1);           % Tonic spikes STA
STAkernel = mean(stimStack,1);
filtStim = cconv(s,STAkernel,length(s));
histFig = figure('Name','Prior and conditional');
%ax(1) = subplot(1,2,1);ps = 
hs = histogram(filtStim,64,'Normalization','probability','DisplayStyle','stairs');
hold on;
hss = histogram(filtStim(stimIdxs),64,'Normalization','probability',...
    'DisplayStyle','stairs');legend({'p(s)','p(s|sp)'})
title('Prior probability OR stimuli probability');xlabel('X (mV)')
ylabel('p(s)');
%ax(2) = subplot(1,2,2);pss = histogram(ax(2),filtStim(stimIdxs),64);
%title('Conditional probability p(stimuli|spike)');xlabel('X (mV)')
%ylabel('p(s|sp)')


STSTDb = std(stimStack(1:Nb,:),[],1);
STSTDt = std(stimStack(Nb+1:end,:),[],1);
end

function [pB,pT,pS] = estimateStimuliPDF(burstStim,tonicStim,stim)
M = 3;
err = 1e-3;
xLow = min(stim);
xHigh = max(stim);
dataRange = [xLow,xHigh];
prB = emforgmm(burstStim,M,err,0);        % Stimuli PDF for busrts
prT = emforgmm(tonicStim,M,err,0);  % Stimuli PDF for tonic sp
prS = emforgmm(stim,M,err,0);        % Stimuli PDF for non-sp
pB = gmmpdf(prB,dataRange);
pT = gmmpdf(prT,dataRange);
pS = gmmpdf(prS,dataRange);
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
function [STAb, STAt] = getSTA(sp,s,maxISI,win,fs)
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
% [stimStack, stimNotSpike] = getEventTriggeredStimulus(winel,Ns,spT,s);
[stimStack, ~] = getEventTriggeredStimulus(winel,Ns,spT,s);

%% Gaussian Mixture Model estimation
% [pB, pT, pS]=...       Bursts            Tonic Spikes          Non-spiking
%     estimateStimuliPDF(stimStack(1:Nb,:),stimStack(Nb+1:end,:),stimNotSpike);
% KL_SB = pS.comparePDF(pB);
% KL_ST = pS.comparePDF(pT);
% Use the histogram values instead of estimating a pdf!!

STAb = mean(stimStack(1:Nb,:),1);               % Burst spikes STA
STAt = mean(stimStack(Nb+1:end,:),1);           % Tonic spikes STA
STAkernel = mean(stimStack,1);
% filtStim = cconv(s,STAkernel,length(s));
% histFig = figure('Name','Prior and conditional');
% hs = histogram(filtStim,64,'Normalization','probability','DisplayStyle','stairs');
% hold on;
% hss = histogram(filtStim(stimIdxs),64,'Normalization','probability',...
%     'DisplayStyle','stairs');legend({'p(s)','p(s|sp)'})
% title('Prior probability OR stimuli probability');xlabel('X (mV)')
% ylabel('p(s)');

STSTDb = std(stimStack(1:Nb,:),[],1);
STSTDt = std(stimStack(Nb+1:end,:),[],1);
end

function [Sstack, NSstack] = getEventTriggeredStimulus(prevSamples,...
    totalSpikes,spikeTimes,stimulus)

Sstack = zeros(totalSpikes,prevSamples);
NSstack = zeros(1,totalSpikes*prevSamples);
suma = 1;
for csp = 1:totalSpikes
    if spikeTimes(csp)-prevSamples+1 >= 1
    Sstack(csp,:) = stimulus(spikeTimes(csp)-prevSamples+1:spikeTimes(csp));
    else, continue;
    end
    if csp == 1
        aux = stimulus(1:spikeTimes(csp)-prevSamples);
    else
        aux = stimulus(spikeTimes(csp-1)+1:spikeTimes(csp)-prevSamples);
    end
    lngth = length(aux);
    NSstack(suma:suma+lngth-1) = aux;
    suma = suma + lngth;
end

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
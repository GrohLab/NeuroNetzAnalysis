timeLapse = [1,5];
fs = 2e4;
load('\\desktop-vui5m85\users\becki\Dropbox\Anesthetized_CT\Data\PairedL5BPOm\130626\p3\freeanalysis.mat','spikeFind*','EEG','*T')
% Creation of the discrete and the continuous stacks:
%                                        First spikeing times  trigger  timeW     FS  FSt  N spiking times          continuous data. 
[discreteStack, contStack] = getStacks(spikeFindingData.spikes,upT,'on',timeLapse,fs,fs,spikeFindingData2.spikes,EEG.data);
%                                                               indx                   cell of indexes | INDX,  signal1, signal2,...
[Ne, Nt, Na] = size(discreteStack);
kIdx = false(1,Na);
binSz = 10; % milliseconds or seconds
ERASE_kIDX = false;
[PSTH, trig, sweeps, tx_psth] = getPSTH(discreteStack, timeLapse, kIdx, binSz, fs);
[triggeredAverageSignals, signalVariation, tx_tav] = getTriggeredAverage(contStack, kIdx, timeLapse);
[relativeSpikeTimes, tx_raster] = getRasterFromStack(discreteStack, kIdx, timeLapse, fs, ERASE_kIDX);

figure;plot(tx_psth,PSTH(1,:),tx_psth,PSTH(2,:))
plotRaster(relativeSpikeTimes, timeLapse, fs, '130626_p3');
plotRasterAsText(relativeSpikeTimes, timeLapse, fs, '130626_p3')
%% Select data file
defaultPath = 'E:\Data\Jesus Emilio\jesus'; % Path where you would like to open a windows explorer window to select your data file
[fileName, filePath] =...
    uigetfile({'*.smr','*.smrx'},'Select a data file',defaultPath);
if fileName == 0
    fprintf(1,'Cancel button pressed. See you next time!\n')
    return
end
%%
timeLapse = [0.05,0.15]; % Time window surrounding the trigger [time before, time after] in seconds
fs = 2e4; % Sampling frequency normally 20 kHz

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
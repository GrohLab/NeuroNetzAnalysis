k = 1e3;m=1e-3;


spkSubs = cellfun(@round,cellfun(@times,sortedData(:,2),...
    repmat({fs},numel(sortedData(:,2)),1),...
    'UniformOutput',false),'UniformOutput',false);
spkLog = cell(size(sortedData,1),1);
for cspk = 1:length(sortedData)
    spkLog(cspk) = {DiscreteWaveform.subs2idx(spkSubs{cspk},Ns)};
end
Triggers = structfun(@double,Triggers,'UniformOutput',false);
%%
timeLapse = [210,150]*m;
binSz = 2*m;
Ncon = numel(Conditions);
for ccon = 1:Ncon
    % Creation of the discrete and the continuous stacks:
    %                                        First spikeing times  trigger  timeW     FS  FSt  N spiking times          continuous data.
    [dst, cst] = getStacks(spkLog{1},Conditions{ccon}.Triggers,'on',timeLapse,fs,fs,...
        spkLog(2:end), struct2cell(Triggers));
    %                                                               indx                   cell of indexes | INDX,  signal1, signal2,...
    [Ne, Nt, Na] = size(dst);
    kIdx = false(1,Na);
    koIdx = true(Ne-2,1);
    tIdx = ~koIdx;
    tIdx(1) = true;
    ERASE_kIDX = false;
    [PSTH, trig, sweeps, tx_psth] = getPSTH(dst, timeLapse, kIdx, binSz, fs);
    [triggeredAverageSignals, signalVariation, tx_tav] = getTriggeredAverage(cst, kIdx, timeLapse);
    [relativeSpikeTimes, tx_raster] = getRasterFromStack(dst, kIdx, koIdx, timeLapse, fs, ERASE_kIDX);
    
    figName = ['M1_',Conditions{ccon}.name,'_S1'];
    plotPSTH(trig,PSTH,sweeps,binSz,timeLapse,figName,...
        cat(1,{Conditions{ccon}.name},sortedData(:,1)),koIdx,tIdx,fs)
    plotRaster(relativeSpikeTimes, timeLapse, fs, figName, sortedData(:,1));
%     plotRasterAsText(relativeSpikeTimes, timeLapse, fs, figName)
end
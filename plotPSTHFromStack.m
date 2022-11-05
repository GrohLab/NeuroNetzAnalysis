function h = plotPSTHFromStack(discStack,trialFlags, timeLapse, binSz, fs,...
    consCondNames)
[Ne, Nt, Na] = size(discStack);
tx = (0:Nt-1)/fs + timeLapse(1);
Na = sum(trialFlags,1);
Np = diff(timeLapse)/binSz;
txPSTH = (0:Np-1)*binSz + timeLapse(1);
h = 
for ccon = 1:size(trialFlags,2)
    PSTH = sum(discStack(:,:,trialFlags(:,ccon)),3);
    tmVals = arrayfun(@(x,y) repmat(x,y,1), tx, PSTH,'UniformOutput',0);
    tmVals = cat(1,tmVals{:});
    h = histogram(tmVals,'BinWidth',binSz,'BinLimits',timeLapse,'Visible','off'); 
    bar(txPSTH, h.Values/(binSz*Na(ccon)),'EdgeColor','none',...
        'FaceAlpha',0.4,'DisplayName',consCondNames{ccon})
    h.delete;
    if ccon == 1
        hold on
    end
end
function plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,IDe,expName)
PSTHn = PSTH./max(PSTH,[],2);
fig = figure('Name',expName,'Color',[1,1,1]);
ax1 = subplot(4,1,1:3,'Parent',fig);
psthTX = linspace(-timeLapse(1),timeLapse(2),size(PSTHn,2));
imagesc(ax1,'XData',psthTX,'CData',PSTHn);
ax1.YLim = [0.5,size(PSTH,1)+0.5];
ax1.XLim = [-timeLapse(1)-binSz/2, timeLapse(2)+binSz/2];
trigTX = linspace(-timeLapse(1),timeLapse(2),size(trig,2));
ax2 = subplot(4,1,4,'Parent',fig);
popPSTH = sum(PSTHn,1)/size(PSTH,1);
plot(ax2,psthTX,popPSTH,'Color',[0.8,0.8,0.8])
end
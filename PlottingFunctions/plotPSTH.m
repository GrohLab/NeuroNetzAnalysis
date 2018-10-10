function ax =...
    plotPSTH(trig, PSTH, sweeps, binSz, timeLapse, expName, IDe, koIdx, tIdx, fs, ax)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
koIdx = koIdx(~tIdx);
trigID = IDe{tIdx};
IDe(tIdx) = [];
if ~exist('ax','var') || isempty(ax)
    fig = figure('Name',['PSTH for ',expName],'Color',[1,1,1]);
    ax = gca;
else
    fig = get(ax,'Parent');
end
tx_PSTH = linspace(-timeLapse(1),timeLapse(2),size(PSTH,2));
tx_trig = linspace(-timeLapse(1),timeLapse(2),length(trig));
set(fig,'defaultAxesColorOrder',[0,0,0;0.2,0.2,0.2])
plot(ax,tx_trig,trig,'DisplayName',trigID,'Color',[37, 154, 3]/255);
ylabel(ax,'Stimulus probability')
yyaxis(ax,'right');
clrMap = jet(size(PSTH,1)-1);
bar(ax,tx_PSTH,PSTH(1,:)/(sweeps*binSz),1,...
    'EdgeColor','none','FaceColor',[0.2,0.2,0.2],...
    'FaceAlpha',0.3,'DisplayName','Neuron 1');
ylabel(ax,'Firing rate [Hz]')
hold(ax,'on')
idIdx = find(koIdx)';
binEl = fs * binSz;
yyaxis(ax,'left')
xlabel(ax,'Time [s]')
for cp = idIdx
    if max(PSTH(cp+1,:)) > sweeps
        plot(ax,tx_PSTH,(PSTH(cp+1,:))/(binEl*sweeps),...
            'Color',clrMap(cp,:),'DisplayName',IDe{cp});
    else
        plot(ax,tx_PSTH,PSTH(cp+1,:)./(binSz*sweeps),...
            'Color',clrMap(cp,:),'DisplayName',IDe{cp});
    end
end
title(ax,[expName,' ',num2str(sweeps),' trials'],'Interpreter','none')
legend(ax,'show')
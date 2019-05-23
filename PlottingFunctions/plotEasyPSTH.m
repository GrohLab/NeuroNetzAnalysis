function fig =...
    plotEasyPSTH(trig, PSTH, sweeps, binSz, timeLapse, fs)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

fig = figure('Name','PSTH','Color',[1,1,1]);
ax = axes('Parent',fig);

tx_PSTH = linspace(-timeLapse(1),timeLapse(2),size(PSTH,2));
tx_trig = linspace(-timeLapse(1),timeLapse(2),length(trig));
plot(ax,tx_trig, trig, 'Color', [37, 154, 3]/255);
ylabel(ax,'Stimulus probability')
yyaxis(ax,'right');
clrMap = jet(size(PSTH,1)-1);
bar(ax,tx_PSTH,PSTH(1,:)/(sweeps*binSz),1,...
    'EdgeColor','none','FaceColor',[0.1,0.1,0.1],...
    'FaceAlpha',0.7);
ylabel(ax,'Firing rate [Hz]')
hold(ax,'on')
ax.Title.Interpreter = 'none';
yyaxis(ax,'left')
xlabel(ax,sprintf('Time_{%.2f} [s]',binSz))
binEl = fs * binSz;
for cp = 1:size(PSTH,1)-1
    if max(PSTH(cp+1,:)) > sweeps
        plot(ax,tx_PSTH,(PSTH(cp+1,:))/(binEl*sweeps),...
            'Color',clrMap(cp,:));
    else
        plot(ax,tx_PSTH,PSTH(cp+1,:)./(binSz*sweeps),...
            'Color',clrMap(cp,:));
    end
end
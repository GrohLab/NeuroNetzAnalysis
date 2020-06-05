function ax =...
    plotPSTH(trig, PSTH, sweeps, binSz, timeLapse, expName, IDe, koIdx, tIdx, fs, ax)
%PLOTPSTH takes the output from the function getPSTH and plots the
%normalized counts of the discrete stack events as either bar plots or as
%lines, depending on their nature. If the event is pulsed, then the plot
%will be a line. Otherwise, it would be a bar i.e. spike times.
%   [axes] = plotPSTH(trig, PSTH, sweeps, binSz, timeLapse, expName, IDe,
%       koIdx, tIdx, fs, ax)
%   
%   INPUTS:
%       trig - Trigger signal. This is an output from the getPSTH function
%       PSTH - Cell array containing the relative times of the events
%       sweeps - Number of trials in the experiment
%       binSz - Bin size (in seconds)
%       timeLapse - Viewing window for the PSTH (in seconds)
%       expName - Name of the experiment
%       IDe - Identification of each event in the PSTH plot
%       koIdx - 
%       tIdx - Logical array indicating the trigger signal (true)
%       fs - Sampling frequency
%       [ax] - Optional axis variable to plot the PSTH in a predefined
%       figure
%   
%   OUTPUTS:
%       ax - axis variable in which the PSTH plot is embedded.
% Emilio Isaias-Camacho @GrohLab 2018

koIdx(tIdx) = false;
trigID = IDe{tIdx};
IDe(tIdx) = [];
if ~exist('ax','var') || isempty(ax)
    fig = figure('Name',['PSTH for ',expName],'Color',[1,1,1]);
    ax = gca;
else
    fig = get(ax,'Parent');
end
tx_PSTH = linspace(timeLapse(1),timeLapse(2),size(PSTH,2));
tx_trig = linspace(timeLapse(1),timeLapse(2),length(trig));
set(fig,'defaultAxesColorOrder',[0,0,0;0.2,0.2,0.2])
plot(ax,tx_trig,trig,'DisplayName',trigID,'Color',[37, 154, 3]/255);
ylabel(ax,'Stimulus probability')
yyaxis(ax,'right');
clrMap = jet(size(PSTH,1));
bar(ax,tx_PSTH,PSTH(1,:)/(sweeps*binSz),1,...
    'EdgeColor','none','FaceColor',[0.2,0.2,0.2],...
    'FaceAlpha',0.3,'DisplayName',IDe{1});
ylabel(ax,'Firing rate [Hz]')
hold(ax,'on')
idIdx = find(koIdx);
idIdx = reshape(idIdx,1,numel(idIdx));
binEl = fs * binSz;
yyaxis(ax,'left')
xlabel(ax,sprintf('Time_{%.3f} [s]',binSz))
for cp = idIdx
    if max(PSTH(cp,:)) > sweeps
        yyaxis(ax,'left')
        plot(ax,tx_PSTH,(PSTH(cp,:))/(binEl*sweeps),...
            'Color',clrMap(cp,:),'DisplayName',IDe{cp});
    else
        yyaxis(ax,'right')
        bar(ax,tx_PSTH,PSTH(cp,:)./(binSz*sweeps),...
            'EdgeColor','none',...
            'FaceColor',clrMap(cp,:),'DisplayName',IDe{cp},...
            'FaceAlpha',0.2);
    end
end
title(ax,[expName,' ',num2str(sweeps),' trials'],'Interpreter','none')
yyaxis('left')
axis(ax,[timeLapse(1),timeLapse(2),0,1.1])
legend(ax,'show','Location','best')
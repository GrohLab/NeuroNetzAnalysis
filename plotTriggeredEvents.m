function [PSTH, figID] =...
    plotTriggeredEvents(PSTHstack, LFPstack, Wstack, timeLapse, cellType, binSz, fs)
% plotTriggeredEvents returns a peri-stimulus triggered histogram given a
% triggered aligned stack in the form MxNxT, where M is the number of
% channels to align, N is the number of samples (time) and T is the number
% of aligning triggers found in the channel. The kicking out process needs
% to be implemented in a dynamical form. The events after the (first) spike
% channel are a guide to kick row of the spikes out. For now, we will kick
% out all of the rows which contain a true in the observed time.

% Computing the size of the PSTH stack
if binSz >= 1
    warning('Assuming binsize given in millisecons. Changing to seconds')
    binSz = binSz*1e-3;
end

[Ne, ~, Na] = size(PSTHstack);
fsLFP = 1e3;
%% Cleaning the stack
% Kicking out everything that contains a true in the channel. No light, no
% puff, no touch. Nothing.
if size(PSTHstack,1) > 2
    kickOutIdx = squeeze(sum(PSTHstack(3:end,:,:),2)) > 0;
    kickAlignmentIDx = sum(kickOutIdx,1) > 0;
else
    kickAlignmentIDx = true(1,Na);
end
disp(['Found extra events: ',num2str(Ne-2)])
disp([num2str(sum(kickAlignmentIDx)), ' omitted alignment points.'])
[PSTH, trig, sweeps] = getPSTH(PSTHstack, kickAlignmentIDx, binSz, fs);
[LFPmean, LFPstd] = getTriggeredAverage(LFPstack, kickAlignmentIDx);
[Wmean, Wstd] = getTriggeredAverage(Wstack, kickAlignmentIDx);
switch cellType
    case 'POm'
        % Black for POm
        colr = 'k';
    case 'VPM'
        % Gray for VPM
        colr = 0.6 * ones(1,3);
    otherwise
        % Orange for other.
        colr = [230, 138, 0]/255;
end
%% Plot the information extracted from the stack.
% Create figure for the PSTH. The color is determined by the cell type
% RASTER PLOT
figID = figure('Name',[char(cellType), ' Raster, PSTH and average traces'],...
    'Color',[1,1,1]);
ax(1) = subplot(5,1,[2,3],'Color','none');
spksStack = squeeze(PSTHstack(2,:,~kickAlignmentIDx));
xx = -timeLapse(1):1/fs:timeLapse(2);
plot(ax(1),[-timeLapse(1),timeLapse(2)],[1,Na-sum(~kickAlignmentIDx)],...
    'LineStyle','none','Marker','none');
for cl = size(spksStack,2):-1:1
    xspks = xx(spksStack(:,cl));
    for cs = 1:numel(xspks)
        text(xspks(cs),cl,'|','HorizontalAlignment','center')
    end
end
ax(1).XColor = 'none';
ylabel('Trials');box off
T = xx(find(PSTHstack(1,:,1)==1,1,'first'));set(gca,'XTickLabel',[])
hold on;plot([T,T],[0,Na-sum(~kickAlignmentIDx)],'Color',[37, 154, 3]/255)
axis([-timeLapse(1),timeLapse(2),0,sum(~kickAlignmentIDx)+1])
% PSTH PLOT
ax(2) = subplot(5,1,4,'Color','none');
bar(ax(2),-timeLapse(1):binSz:timeLapse(2),PSTH/(binSz*sweeps),1,...
    'EdgeColor','none','FaceColor',colr);
hold on;ylabel('Frequency [Hz]');box off;xlabel(['Time_{',num2str(binSz*1e3),' ms} [s]'])
binEls = round(binSz * fs);
yyaxis(ax(2),'right')
area(ax(2),-timeLapse(1):binSz:timeLapse(2),(trig/(binEls*sweeps)),...
    'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[37, 154, 3]/255)
ylabel('Probability of whisking');set(gca,'YColor',[37, 154, 3]/255)
% AVERAGE LFP
axLFP(1) = subplot(5,1,1);
plotAverageTrace(axLFP(1),timeLapse,LFPmean,LFPstd,'LFP',fsLFP)

% AVERAGE WHISKER
axLFP(2) = subplot(5,1,5);
plotAverageTrace(axLFP(2),timeLapse,Wmean,Wstd,'Whisker',fsLFP)
linkaxes([ax,axLFP],'x')
end

function plotAverageTrace(ax,timeLapse,y,sigma,signalName,fs)
txLFP = -timeLapse(1):1/fs:timeLapse(2);
fill(ax,[txLFP,fliplr(txLFP)],[y'-sigma',fliplr(y'+sigma')],...
    [0.9,0.9,0.9],'LineStyle','none');hold on
plot(ax,txLFP,y,'LineWidth',3,'Color',[99, 115, 131]/255)
line(ax,[0,0],[min(y-sigma),max(y+sigma)],'Color',[37, 154, 3]/255)
text(-timeLapse(1),mean(y),signalName,'FontWeight','bold',...
    'HorizontalAlignment','right')
ax.Color = 'none';
ax.XColor = 'none';
ax.YColor = 'none';
end
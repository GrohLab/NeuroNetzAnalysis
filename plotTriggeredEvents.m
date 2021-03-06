function [PSTH, figID, kickAlignmentIDx] =...
    plotTriggeredEvents(PSTHstack, LFPstack, Wstack, timeLapse, cellType, binSz, fs, fsLFP)
% PLOTTRIGGEREDEVENTS returns a peri-stimulus triggered histogram given a
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

[Ne, Nt, Na] = size(PSTHstack);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEWARE OF THE LFP SAMPLING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FREQUENCY
% fsLFP = 1e3;
%% Cleaning the stack
% Kicking out everything that contains a true in the channel. No light, no
% puff, no touch. Nothing. Inverted variable: kick means to keep.
if size(PSTHstack,1) > 2
    kickOutIdx = squeeze(sum(PSTHstack(3:end,:,:),2)) > 0;
    kickAlignmentIDx = sum(kickOutIdx,1) > 0;
else
    kickAlignmentIDx = false(1,Na);
end
disp(['Found extra events: ',num2str(Ne-2)])
disp([num2str(sum(kickAlignmentIDx)), ' omitted alignment points.'])
[PSTH, trig, sweeps] = getPSTH(PSTHstack, timeLapse, kickAlignmentIDx, binSz, fs);
[LFPmean, LFPstd] = getTriggeredAverage(LFPstack, kickAlignmentIDx, timeLapse);
[Wmean, Wstd] = getTriggeredAverage(Wstack, kickAlignmentIDx, timeLapse);
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
figID(1) = figure('Name',[char(cellType), ' Raster, PSTH and average traces'],...
    'Color',[1,1,1]);
ax(1) = subplot(5,1,[2,3],'Color','none');
spksStack = squeeze(PSTHstack(2,:,~kickAlignmentIDx));
tx = 0:1/fs:(length(trig)-1)/fs;
tx = tx - timeLapse(1);
plot(ax(1),[-timeLapse(1),timeLapse(2)],[1,Na-sum(~kickAlignmentIDx)],...
    'LineStyle','none','Marker','none');
% Plot spikes as text.
for cl = size(spksStack,2):-1:1
    xspks = tx(spksStack(:,cl));
    for cs = 1:numel(xspks)
        text(xspks(cs),cl,'|','HorizontalAlignment','center')
    end
end
ax(1).XColor = 'none';
ylabel(['Trials (',num2str(sum(~kickAlignmentIDx)),'/',num2str(Na),')']);
box off;T = tx(round((Nt*timeLapse(1) - timeLapse(2))/sum(timeLapse))+1);
% T = xx(find(PSTHstack(1,:,1)==1,1,'first'));set(gca,'XTickLabel',[])
set(gca,'XTickLabel',[]);hold on;
plot([T,T],[0,sum(~kickAlignmentIDx)],'Color',[37, 154, 3]/255)
axis([-timeLapse(1),timeLapse(2),0,sum(~kickAlignmentIDx)+1])
% PSTH PLOT
ax(2) = subplot(5,1,4,'Color','none');
bar(ax(2),-timeLapse(1):binSz:timeLapse(2),PSTH(1,:)/(binSz*sweeps),1,...
    'EdgeColor','none','FaceColor',colr);
hold on;ylabel('Frequency [Hz]');box off;xlabel(['Time_{',num2str(binSz*1e3),' ms} [s]'])
% binEls = round(binSz * fs);
yyaxis(ax(2),'right')
% area(ax(2),-timeLapse(1):binSz:timeLapse(2),(trig/(binEls*sweeps)),...
%     'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[37, 154, 3]/255)

%area(ax(2),-timeLapse(1):1/fs:timeLapse(2),trig,...
%    'EdgeColor','none','FaceAlpha',0.3,'FaceColor',[37, 154, 3]/255)
area(ax(2),tx,trig,...
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
if ~isempty(y)
    y = double(y);
    tx = 0:1/fs:(length(y)-1)/fs;
    tx = tx - timeLapse(1);
    fill(ax,[tx,fliplr(tx)],[y-sigma,fliplr(y+sigma)],...
        [0.9,0.9,0.9],'LineStyle','none');hold on
    plot(ax,tx,y,'LineWidth',3,'Color',[99, 115, 131]/255)
    line(ax,[0,0],[min(y-sigma),max(y+sigma)],'Color',[37, 154, 3]/255)
    text(tx(1),mean(y),signalName,'FontWeight','bold',...
        'HorizontalAlignment','right')
end
ax.Color = 'none';
ax.XColor = 'none';
ax.YColor = 'none';
end
function [fig] =...
    plotRaster(relativeSpikeTimes, timeLapse, fs, figTitle, IDe, fig)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Ne is the number of 'events', namely signals that are considered in the
% stack, i.e. neurons
% Na is the number of alignment points in this stack
[Ne, Na] = size(relativeSpikeTimes);
% Time axis
idxOffset = timeLapse * fs;
Nt = sum(idxOffset) + 1;
tx = 0:1/fs:(Nt-1)/fs;
tx = tx - timeLapse(1);
if ~exist('figTitle','var')
    figTitle = [];
end
plotTitle = ['Raster Plot ', figTitle, ' ',num2str(Na),' trials'];
AX_FLAG = true;
if ~exist('fig','var') || isempty(fig)
    fig = figure('Name',plotTitle,'Color',[1,1,1]);
    AX_FLAG = false;
end
cmap = [0.01,0.01,0.01;jet(Ne - 1)];
FIRST_FLAG = true;
yTickLabel = cell(1,Ne);
for cse = 1:Ne
    % For each spike train
    yTickLabel(cse) = IDe(cse);
    for cap = Na:-1:1
        % For each alignment point
        if ~isempty(relativeSpikeTimes{cse,cap})
            xspks = tx(relativeSpikeTimes{cse,cap} + idxOffset(1) + 1);
            lvl = (cse - 1)*Na + cap;
            if FIRST_FLAG
                if AX_FLAG
                    hold on
                else
                    hold on
                end
                FIRST_FLAG = false;
            end
            if AX_FLAG
                plot(xspks,lvl*ones(1,numel(xspks)),...
                    'LineStyle','none','Marker','.',...
                    'MarkerFaceColor',cmap(cse,:),'MarkerSize',2);
            else
                plot(xspks,lvl*ones(1,numel(xspks)),...
                    'LineStyle','none','Marker','.',...
                    'MarkerFaceColor',cmap(cse,:),'MarkerSize',2,...
                    'Color',cmap(cse,:))
            end
            
        end
    end
end

axis([-timeLapse(1),timeLapse(2),1,Na*Ne])
if AX_FLAG
    set(fig,'Box','off','YTick',(0:Ne)*Na + Na/2,'YTickLabel',yTickLabel,...
        'YTickLabelRotation',90)
    xlabel(fig,'Time [s]')
    title(fig,[figTitle, ' ',num2str(Na),' trials'],'Interpreter','none')
else
    fig = gca;
    set(fig,'Box','off','YTick',(0:Ne)*Na + Na/2,'YTickLabel',yTickLabel,...
        'YTickLabelRotation',90)
    xlabel('Time [s]')
    title([figTitle, ' ',num2str(Na),' trials'],'Interpreter','none')
end

end


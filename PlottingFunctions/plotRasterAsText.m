function [ax] =...
    plotRasterAsText(relativeSpikeTimes, timeLapse, fs, figTitle, ax)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Ne is the number of 'events', namely signals that are considered in the
% stack
% Nt is the time axis for the chosen window
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
plotTitle = ['Raster Plot ', figTitle];
if ~exist('ax','var') || isempty(ax)
    figure('Name',plotTitle,'Color',[1,1,1]);plot(NaN,NaN);
    ax = gca;
end
cmap = colormap(jet(Ne));
yTickLabel = cell(1,Ne);
for cse = 1:Ne
    % For each spike train
    yTickLabel(cse) = {['Neuron ',num2str(cse)]};
    for cap = Na:-1:1
        % For each alignment point
        if ~isempty(relativeSpikeTimes{cse,cap})
            xspks = relativeSpikeTimes{cse,cap};
            lvl = (cse - 1)*Na + cap;
            for cs = 1:numel(xspks)
                % For each spike
                text(ax, xspks(cs),lvl,'.',...
                    'HorizontalAlignment','center','Color', cmap(cse,:))
            end
        end
    end
end

axis([-timeLapse(1),timeLapse(2),1,Na*Ne])
ax = gca;set(ax,'Box','off','YTick',(0:Ne)*Na + Na/2,'YTickLabel',yTickLabel,...
    'YTickLabelRotation',90)
xlabel(ax,'Time [s]')
title(ax,figTitle,'Interpreter','none')


end


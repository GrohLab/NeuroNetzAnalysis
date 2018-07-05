function [ax] =...
    plotRaster(relativeSpikeTimes, timeLapse, fs, figTitle, ax)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% Ne is the number of 'events', namely signals that are considered in the
% stack, i.e. neurons
% Na is the number of alignment points in this stack
[Ne, Na] = size(relativeSpikeTimes);
% Time axis
tx = 0:1/fs:(Nt-1)/fs;
tx = tx - timeLapse(1);

AX_FLAG = true;
if ~exist('ax','var') || isempty(ax)
    figure('Name',['Raster Plot ', figTitle],'Color',[1,1,1]);
    AX_FLAG = false;
end
cmap = colormap(jet(Ne-1));
FIRST_FLAG = true;
idxOffset = timeLapse(1) * fs;
for cse = 2:Ne
    % For each spike train
    for cap = Na:-1:1
        % For each alignment point
        xspks = tx(relativeSpikeTimes{cse,cap} + idxOffset);
        lvl = (cse - 2)*Na + cap;
        if AX_FLAG
            plot(ax,xspks,lvl*ones(1,numel(xspks)),...
                'LineStyle','none','Marker','.',...
                'MarkerFaceColor',cmap(cse-1,:),'MarkerSize',2)
        else
            plot(xspks,lvl*ones(1,numel(xspks)),...
                'LineStyle','none','Marker','.',...
                'MarkerFaceColor',cmap(cse-1,:),'MarkerSize',2,...
                'Color',cmap(cse-1,:))
        end
        if FIRST_FLAG
            hold on
            FIRST_FLAG = false;
        end
    end
end
axis([-timeLapse(1),timeLapse(2),1,Na*(Ne-1)])
ax = gca;
end


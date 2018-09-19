function [outputArg1,outputArg2] =...
    plotExampleTraces(discreteTraces,IDs,plotTitle,fs)
%PLOTEXAMPLETRACES Precisely plots the spikes and the other considered
%events in time to have a first data observation.
%   Detailed explanation goes here
N = size(discreteTraces,2);
tx = 0:1/fs:(N-1)/fs;
fig = figure('Name','Example traces','Color',[1,1,1]);
axes
ax = get(fig,'CurrentAxes');
lbls  = false(size(discreteTraces,1),1);
for ct = size(discreteTraces,1):-1:1
    pFlag = plot(ax,tx(discreteTraces(ct,:)),...
        ct*ones(ct,sum(discreteTraces(ct,:))),...
        'LineStyle','none','Marker','.');
    if ~isempty(pFlag)
        lbls(ct) = true;
    end
    if ct == size(discreteTraces,1)
        hold on
    end
end
lvls = find(lbls);
set(ax,'YTick',lvls,'YTickLabel',IDs(lbls))
xlabel(ax,'Time [s]');title(ax,plotTitle,'Interpreter','none')
axis([0 tx(end) lvls(1)-0.5 lvls(end)+0.5])
end


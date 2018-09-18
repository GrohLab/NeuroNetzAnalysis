function [outputArg1,outputArg2] = plotExampleTraces(discreteTraces,IDs,fs)
%PLOTEXAMPLETRACES Precisely plots the spikes and the other considered
%events in time to have a first data observation.
%   Detailed explanation goes here
N = size(discreteTraces,2);
tx = 0:1/fs:(N-1)/fs;
fig = figure('Name','Example traces','Color',[1,1,1]);
axes
ax = get(fig,'CurrentAxes');
for ct = size(discreteTraces,1):-1:1
    plot(ax,tx(discreteTraces(ct,:)),ct*ones(ct,sum(discreteTraces(ct,:))),'LineStyle','none','Marker','.')
    if ct == size(discreteTraces,1)
        hold on
    end
end
set(ax,'YTickLabel',IDs)
end


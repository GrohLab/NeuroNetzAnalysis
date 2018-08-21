function ax =...
    plotTriggeredAverage(avSignals, avSTD, tx, IDs,expName, ax)
%UNTITLED5 Summary of this function goes here  triggeredAverageSignals,signalVariation,tx,IDs
%   Detailed explanation goes here
if ~exist('ax','var') || isempty(ax)
    figure('Name',['Triggered average for ',expName],'Color',[1,1,1])
    ax = gca;
    close gcf
end
Ne = size(avSignals,1);
clrMap = jet(Ne);
FF = true;
if Ne <= 2
    figure('Name',['Triggered average for ',expName],'Color',[1,1,1])
    ax = gca;plot(ax,tx,tx>=0,'Color',[37, 154, 3]/255,'DisplayName','Trigger')
    yyaxis('right');plot(ax,tx,avSignals(1,:),'DisplayName',IDs{1});
    legend('show')
    figure('Name',['Triggered average for ',expName],'Color',[1,1,1])
    ax = gca;plot(ax,tx,tx>=0,'Color',[37, 154, 3]/255,'DisplayName','Trigger')
    yyaxis('right');plot(ax,tx,avSignals(2,:),'DisplayName',IDs{2})
    legend('show')
else
    for cs = 1:size(avSignals,1)
        plot(ax,tx,avSignals(cs,:),'Color',clrMap(cs,:),'DisplayName',IDs{cs})
        if FF
            hold on;FF = false;
        end
    end
end
legend('show')
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
    yyaxis('right');pl1= plot(ax,tx,avSignals(1,:),'DisplayName',IDs{1});
    hold on;plot(ax,tx,avSTD(1,:)+avSignals(1,:),'LineStyle','--')
    plot(ax,tx,avSignals(1,:)-avSTD(1,:),'LineStyle','--')
    legend(ax,pl1,IDs{1})
    figure('Name',['Triggered average for ',expName],'Color',[1,1,1])
    ax = gca;plot(ax,tx,tx>=0,'Color',[37, 154, 3]/255,'DisplayName','Trigger')
    yyaxis('right');pl2 = plot(ax,tx,avSignals(2,:),'DisplayName',IDs{2});
    hold on;plot(ax,tx,avSTD(2,:)+avSignals(2,:),'LineStyle','--')
    plot(ax,tx,avSignals(2,:)-avSTD(2,:),'LineStyle','--')
    legend(ax,pl2,IDs{2})
else
    for cs = 1:size(avSignals,1)
        pl(cs) = plot(ax,tx,avSignals(cs,:),'Color',clrMap(cs,:),'DisplayName',IDs{cs});
        if FF
            hold on;FF = false;
        end
        plot(ax,tx,avSTD(cs,:)+avSignals(cs,:),'Color',clrMap(cs,:))
        plot(ax,tx,avSignals(cs,:)-avSTD(cs,:),'Color',clrMap(cs,:))
    end
    legend(pl,'show')
end

function [sp h bins uptrig_EEG uptrig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,str, EEG,light, plotit, plotRaster);
if ~exist('plotRaster','var')
    plotRaster = true;
end

EEG=EEG-min(EEG);EEG=EEG/max(EEG);
light=light-min(light);light=light/max(light);

f=0;
if size(spikes,2)>size(spikes,1),spikes=spikes';end
sp={};
for i=1:numel(triggers)
    spnorm=spikes-triggers(i);
    sp{i}=spnorm(find(spnorm>=-timeBefore & spnorm<=timeAfter));
end

if plotit
    f=figure;
    
    if plotRaster
        subplot(4,1,1:2)
        %ylim([0 numel(sp)+1])
        for i=1:numel(sp)
            
            plot(sp{i}/ppms,ones(size(sp{i}))*i,'k.')
            %         %ram quick fix to deal with weird raster issue
            %         if ~isempty(sp{i})
            %             r=raster(sp{i}/ppms,i,'k',.9)
            %         end
            hold on
            
            
        end
        %     try
        %         set(r,'ShowBaseLine','off');
        %     catch
        %     end
        ylabel('trial')
        ylim([0 numel(sp)])
        xlim([-timeBefore/ppms  timeAfter/ppms])
    end
    
end
bins=[-timeBefore:binsize:timeAfter]/ppms;
h=hist(cell2mat(sp')/ppms,bins)/numel(sp);
if plotRaster && plotit
    subplot(4,1,3:4)
end
if plotit
    plot([0 0], [0 numel(sp)+1],'r:','linewidth',3)
    title(str)

    bar(bins,h,1,'k');
    xlim([-timeBefore/ppms  timeAfter/ppms])
    hold on
    plot([0 0], get(gca,'ylim'),'r:','linewidth',3)
    
    ylabel('count/trial')
    
end


uptrig_EEG=0;
if ~isempty(EEG)
    uptrig_EEG=TriggeredSegments(EEG,triggers,-timeBefore, timeAfter)';
end


uptrig_light=0;
if ~isempty(light)
    uptrig_light=TriggeredSegments(light,triggers,-timeBefore, timeAfter)';
end



if plotit
    t=[1:size(uptrig_light,2)]/ppms-timeBefore/ppms;
    m=get(gca,'ylim');m=m(end)*.9;
    p=plot(t,mean(uptrig_light)*m,'b','linewidth',.5)
    hold on
    axis tight
    ylim(get(gca,'ylim')*1.5);
    plot([0 0], get(gca,'ylim'),'r:','linewidth',3)
    xlabel('ms relative to trigger')
    set(gcf,'position',[597    49   483   636])
end


if plotit
    m=get(gca,'ylim');m=m(end)*.9;
    p=plot(t,mean(uptrig_EEG)*m,'g','linewidth',2)
    hold on
    axis tight
    ylim(get(gca,'ylim')*1.5);
    plot([0 0], get(gca,'ylim'),'r:','linewidth',3)
    xlabel('ms relative to trigger')
    set(gcf,'position',[597    49   483   636])
end

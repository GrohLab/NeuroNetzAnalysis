%% multiunit recording practice
clear all
cd 'D:\Dropbox\16 Channel Recording may 2018'
clear all
load SpikeTimes_all_channels.mat
load M137_C5_Mech_L6_05mW_Triggersanalysis Conditions Triggers



%%

%population PSTHs
Spikes={};
Names={};

for i=1:size(sortedData,1)
    Spikes{i}=cell2mat(sortedData(i,2));
    Names{i}=sortedData(i,1);
end

mech=Triggers.whisker;
light=Triggers.light;

close all
bads=[1 2 3 14 6]  %1 2 3 14 are light artifacts
noresponse=[15]
bads=[bads noresponse]
goods=1:16;goods=setdiff(goods,bads);

goods=([16 7 8 12 13 10 9 11 4 5 ]); % make raster plot so that top plot is most superficial neuron
Spikes{5}=sort([Spikes{5} Spikes{6}]);
%goods=[9 4 11 5 7 8]
Spikes=Spikes(goods)
Names=Names(goods);
%% looking at collected data
for I=[1 3:4]
    ppms=20;
    spikes=cell2mat(Spikes)*1000*ppms; %spikes back in samples
    name=Names{I};
    timeBefore=200*ppms;timeAfter=5500*ppms;
    plotit=1;binsize=200;
    triggers=Conditions{I}.Triggers;
    [sp h bins trig_mech trig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
end


subplot(2,1)
bar(bins,h,1,'k')
%% looking at individual

close all
for i=1:numel(Spikes)
    if ~ismember(i,bads)
        for I=[3]
            ppms=20;
            spikes=(Spikes{i})*1000*ppms; %spikes back in samples
            name=Names{i};
            timeBefore=1000*ppms;timeAfter=9000*ppms;
            plotit=1;binsize=100*ppms;
            triggers=Conditions{I}.Triggers;
            [sp h bins trig_mech trig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
            title(i)
        end
        
    end
end


%% organize real figure


% population histogram for later plotting, by condition
H=[];
conds={};
count=0;
Trig_mech={};Trig_light={}
for I=[1 3:4]
    count=count+1;
    ppms=20;
    spikes=cell2mat(Spikes)*1000*ppms; %spikes back in samples
    name=Names{I};
    timeBefore=1000*ppms;timeAfter=5500*ppms;
    plotit=0;binsize=50*ppms;
    triggers=Conditions{I}.Triggers;
    
    conds{count}=Conditions{I}.name;
    [sp h bins Trig_mech{count} Trig_light{count} f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
    
    %convert to rate in Hz
    h=h*(1000/binsize*ppms);
    H(count,:)=h;
end
t=[-timeBefore:timeAfter]/ppms;
%
%get individual responses by condition

Sp={};
close all
count=0;
for I=[1 3:4]  %by condition
    count=count+1;
    triggers=Conditions{I}.Triggers;
    SPIKES={};
    for i=1:numel(Spikes)   %for each neuron
        ppms=20;
        spikes=(Spikes{i})*1000*ppms; %spikes back in samples
        name=Names{i};
        plotit=0;
        %use same binning etc as pop hist
        [SPIKES{i} h bins trig_mech trig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
    end
    Sp{count}=SPIKES;
end



%

%get data appropriate for plotting rasters


SPIKESs={};YSs={}; %all data here ,cond x num neuron

for ii=1:numel(Sp); %pick one condition
    sp=Sp{ii};
    shift=0; SPIKES={};YS={};
    for j=1:numel(sp) %over all neurons
        spikes=sp{j};
        ys={};
        for jj=1:numel(spikes) %over all trials
            if ~isempty(spikes{jj})
                ys{jj}=ones(size(spikes{jj}))+shift;
            else
                ys{jj}=[];
            end
            shift=shift+1;%add for each trial, per neuron
        end
        SPIKES{j}=cell2mat(spikes');
        YS{j}=cell2mat(ys');
        
    end
    
    SPIKESs{ii}=SPIKES;
    YSs{ii}=YS;
end
%

%get color for each neuron, just for plotting
colormap jet;
cmap=colormap;
n=floor(size(cmap,1)/numel(SPIKESs{1}));
colors=cmap(1:n:end,:)


%% plot it all
figure
for ii=1:numel(SPIKESs)
    subplot(6,3,[ii ii+3 ii+6])
    
    for j=1:numel(SPIKESs{ii})
        xs=SPIKESs{ii}{j}/ppms;
        ys=YSs{ii}{j}
        plot(xs,ys,'.','color',colors(j,:),'markersize',10)
        hold on
    end
   % ylabel 'trials/neuron'
    box off
    xlim([min(bins) max(bins)])
end


tits={'mechanical', 
    'mechanical + 10 Hz L6', 
    '10 Hz L6 control'}
for ii=1:size(H,1)
    
    subplot(6,3,[ii+12 ii+15])
    bar(bins,H(ii,:),'k')
    xlim([min(bins) max(bins)])
    xlabel ms
   % ylabel 'pooled spike rate'
    title(tits{ii})
    box off
    ylim([0 200])
end

%
subplot(6,3,10)
plot(t,Trig_mech{1}(1,:),'g','linewidth',2)
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])
ylabel Stimulus

subplot(6,3,11)
plot(t,Trig_mech{2}(1,:),'g','linewidth',2);
hold on
plot(t,Trig_light{2}(1,:),'c','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])

subplot(6,3,12)

plot(t,Trig_light{3}(1,:),'c','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])


%% plot some raw data
figure
load M137_C5_Mech_L6_05mW_Triggersanalysis filteredResponse

v=filteredResponse.data;indices=[10000:(1000*ppms*700)]+30*1000*ppms; v=v/max(v);
time=indices/ppms/1000;
light_in=light(indices);light_in=light_in/max(light_in);
mech_in=mech(indices);mech_in=mech_in/max(mech_in);
plot(time,v(indices)/2,'k')
hold on
for i=1:numel(Spikes)
    
    sp=Spikes{i}*1000*ppms;
    sp=sp(ismember(sp,indices));
    plot(sp/ppms/1000, ones(size(sp))+i*.05-.7,'.','color',colors(i,:),'markersize',15)
    hold on
    axis tight
    ylim([-1.4 1])
    
end
plot(time,light_in/4-.8,'c','linewidth',.5)
plot(time,mech_in/4-1.26,'g','linewidth',1.5)
    



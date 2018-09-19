figure;
if ~exist('Conditions','var') || ~exist('Triggers','var')
    try
        load([rootNames,'analysis.mat'],'Conditions','Triggers')
        disp('Loading successful')
    catch
        fileNames = dir('*analysis.mat');
        load([fileNames.name,'analysis.mat'],'Conditions','Triggers')
        disp('Running without context...')
        fprintf('File %s loaded',fileNames.name)
    end
end

Spikes={};
Names={};

for i=1:size(sortedData,1)
    Spikes{i}=cell2mat(sortedData(i,2));
    Names{i}=sortedData(i,1);
end

mech=Triggers.whisker;
light=Triggers.light;
lenSpks = length(Spikes);
try
    tx = 0:1/Fs:(length(mech) - 1)/Fs;
catch
    disp('The variable Fs doesn''t exist in the workspace')
    disp('Displaying indexes instead of time')
    tx = 1:length(mech);
end
ax(1) = subplot(4,1,1);plot(tx,light,tx,mech)
ax(2) = subplot(4,1,2:4);
Ncl = size(sortedData,1);
fs = dataLoader.SamplingFrequency;

lvls = 1;
bads = [];
consIdxs = true(1,lenSpks);
consIdxs(bads) = false;
lbls = cell(1,sum(consIdxs));
for ccl = 1:Ncl
    if consIdxs(ccl)
        lbls(lvls) = {num2str(ccl)};
        plot(tx(round(sortedData{ccl,2}*fs)),...
            lvls*ones(1,numel(sortedData{ccl,2})),'LineStyle','none',...
            'Marker','.')
        if lvls == 1
            hold on
        end
        lvls = lvls + 1;
    end
end
ylabel('Cluster Number');xlabel('Time [s]')
axis([tx(1),tx(2),0,lvls+1])
set(ax(2),'YTick',1:lvls-1,'YTickLabels',lbls)
linkaxes(ax,'x')

Ncond = numel(Conditions);
clrs = [255, 128, 0;...
    0, 255, 255;...
    0, 128, 255;...
    128, 0, 255]/255;
for ccon = 1:Ncond
    Nt = numel(Conditions{ccon}.Triggers);
    plot(ax(2),[tx(Conditions{ccon}.Triggers);tx(Conditions{ccon}.Triggers)],...
        [(1+Ncl)*ones(1,Nt);zeros(1,Nt)],'Color',clrs(ccon,:))
end
hold off
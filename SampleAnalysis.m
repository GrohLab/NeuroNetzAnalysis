%% multiunit recording practice
% cd 'D:\Dropbox\16 Channel Recording may 2018'
% The clearvars command is commented to avoid erasing the information
% extracted from the smrx file through the UMS_life_script.mlx
% clearvars

homedir='C:\Users\NeuroNetz\Documents\Data\Sailaja\2019\Alex\M137_4';
% homedir='E:\Data\';
cd(homedir)
fname = 'M137_C5_Mech_L6 05mW'; 

load([fname,'_all_channels.mat'])
load([fname,'analysis.mat'],'Conditions','Triggers')
try
    fs = Fs;
catch
    disp('Fs is not existing...')
    disp('Defining...')
    Fs  = fs;
end

%% Initialize the variables
% This section loads the cluster spike times into the 'Spikes' cell array.
% If there is a step gone wrong, you can re-initialize re-running this
% section.

%population PSTHs
Spikes={};
Names={};

for i=1:size(sortedData,1)
    Spikes{i}=cell2mat(sortedData(i,2));
    Names(i)=sortedData(i,1);
end

mech=Triggers.whisker;
light=Triggers.light;

close all

%% look at inter-spike intervals as another check for pure light-evoke electrical artifacts
close all
for i=1:numel(Spikes) % insert cluster numbers here
    figure
      isis=[diff(Spikes{i})*1000]; %in ms
      [hisi isi_bins]=hist(isis,[1:10000]);
      subplot(2,1,1)
      plot(log10(isi_bins),hisi,'linewidth',2)
       % title(['Cluster: ',num2str(i)])
      title(['Cluster ID: ',Names{i},' #',num2str(i)],'Interpreter','none')
      subplot(2,1,2)
      abc = string(sortedData(i,1));
      C = strsplit(abc,'_');
      D = C(1,1);
      E = erase(D, 'ch');
      F = C(1,3);
      load([fname,'_',num2str(E),'.mat']);
      plot_waveforms(spikes,str2num(F));
end
        

%% Possible artifacts (is laser evoking a response or an optoelectric artifact?)

%bads=[1 2 3 14 6]  %1 2 3 14 are light artifacts

bads=[12 18 19];

noresponse=[];
bads=[bads noresponse];
%bads=[]; %uncomment this to have bads empty
%assumes that all spike trains are good
goods=1:numel(Spikes);

%removes bad units 
goods=setdiff(goods,bads)

%% Let us see if the laser evokes a neural response.
%use these two lines to look at hand-defined lasers response

%goods=[]; %fill in hand-selected channels here!
%bads=[];

%Spikes=Spikes(goods);
%Names=Names(goods);




%% looking at collected data
fs = Fs;

%close all

for I=[1:4]
    ppms=fs/1000;
    spikes=cell2mat(Spikes)*1000*ppms; %spikes back in samples
    name=Names{I};
    timeBefore=500*ppms;timeAfter=6500*ppms;
    plotit=1;binsize=2.5*ppms;
    triggers=Conditions{I}.Triggers;
    [sp h bins trig_mech trig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
    title(Conditions{I}.name)
end

lenSpks = length(Spikes);
%% Cross-correlation or relationship between clusters
% Creating a square matrix containing the cross correlation normalized
% peaks omitting the 'bads' clusters. The input to the cross correlation
% function should be 'normalized'. This means that the spike traces shall
% be reconstructed from the time stamps and input to the xcorr function.
lenSpks = length(Spikes);
consIdxs = true(1,lenSpks);
% Some clusters were eliminated after an induvidual inspection of their
% PSTHs.
% bads = [];
consIdxs(bads) = false;
crscor = zeros(lenSpks,lenSpks,3);    % Square matrix with Ncl x Ncl
alex = true;
for ccl = 1:lenSpks
    xccl = ccl+1;
    % Cross correlation avoiding the autocorrelations and the 'bads'.
    % AM: Lag-limited cross-correlation and display together with the
    % distance matrix.
    while consIdxs(ccl) && xccl <= lenSpks
        % auxXCOR = xcorr(Spikes{ccl},Spikes{xccl},'coeff');
        if consIdxs(xccl)
            clstrInfo =...
                ['Cluster ',num2str(ccl),' against cluster ',num2str(xccl)];
            disp(clstrInfo)
            % Signal reconstruction
            auxSignal1 = false(1,length(mech));
            auxSignal2 = false(1,length(mech));
            % Spikes assignment
            auxSignal1(round(fs*Spikes{ccl})) = true;
            auxSignal2(round(fs*Spikes{xccl})) = true;
            % Selecting a subsampled version of the signals given a maximum
            % lag.
            MxLag = 25; % Seconds
            rIdx = randi([round(fs*MxLag),length(mech)-round(fs*MxLag)],1);
            disp(['Random time selected: ',num2str(rIdx/fs),' seconds'])
            rWindIdx = rIdx-round(fs*MxLag):rIdx+round(fs*MxLag);
            auxSignal1 = auxSignal1(rWindIdx);
            auxSignal2 = auxSignal2(rWindIdx);
            
            % Distance matrix -- AM: promt the user for keeping the current
            % DM for later processing. (Inter-cluster spike intervals)
            dfMtx = log(distmatrix(Spikes{ccl}',Spikes{xccl}')+1);
            lnIdx = dfMtx < log(16);
            [y,x]=find(lnIdx);
            [mdl,yhat,r2] = fit_poly(x,y,1);
            eqInfo = ['y = ',num2str(mdl(1)),'x ',num2str(mdl(2))];
            display(eqInfo)
            figure('Name',clstrInfo);
            subplot(2,1,1)
            ax4 = subplot(2,1,1)
            colormap(ax4, spring())
            imagesc(dfMtx);hold on;plot(x,yhat,'LineStyle','--',...
                'LineWidth',3,'Color',[0,1,0]);title([eqInfo,' ',num2str(r2)])
            crscor(ccl,xccl,1:2) = mdl;
            crscor(ccl,xccl,3) = r2;
            subplot(2,1,2)
            abc = string(sortedData(xccl,1));
            C = strsplit(abc,'_');
            D = C(1,1);
            E = erase(D, 'ch');
            F = C(1,3);
            load([fname,'_',num2str(E),'.mat']);
            plot_waveforms(spikes,str2num(F));
            % used for saving plot + question box
            cname =...
                   ['Cl_',num2str(ccl),'_vs_Cl_',num2str(xccl)]
            if alex
                answer = questdlg('Do you want to save this comparison?', ['Cluster Comparison ' cname], 'Save', 'No',"Don't ask again", 'No')
                %Handle response to question
                switch answer
                    case 'Save'
                       return
                       %save(cname)
                    case 'No'
                       %remove the return to run it properly, it is
                       %just to break the loop
                       disp(['Not saving: ' cname])
                    case "Don't ask again"
                       disp('test')
                       alex = false;
                end
            else
                disp('no longer asking to save')
            end
        end
        xccl = xccl + 1;
    end
end
% AM: Check out the waveforms of the interesting clusters

%save cross-correlation so that we don't need to redo this constantly
cd(homedir)
% SAVE FOR THE 'NORMAL' PIPELINE
save CrossCoeffData crscor goods bads
%SAVE FOR THE LASER RESPONSE
% save LaserResponse crscor goods bads
%% Merge similar clusters
% The marging packages indicate which clusters should be merged together
% due to their high similarity. The possibilities are that they belong to a
% same unit as busrting spikes, the cell shifted to another channel or any
% other reasonable cause.
mergingPackages = {[16,17]};
Npg = numel(mergingPackages);
mSpikes = cell(1,Npg);
auxSignal = false(1,length(mech));
for cpg = 1:Npg
    for ccl = 1:numel(mergingPackages{cpg})
        auxSignal(round(fs*Spikes{mergingPackages{cpg}(ccl)})) = true;
    end
    mSpikes(cpg) = {find(auxSignal)/fs};
    auxSignal = false(1,length(mech));
    bads = [bads mergingPackages{cpg}(2:end)];
    Spikes(mergingPackages{cpg}(1)) = {single(mSpikes{cpg})};
end
bads = sort(unique(bads))

%save merging info
%SAVE FOR 'NORMAL' PIPELINE
save CrossCoeffData Spikes mergingPackages bads -append
% SAVE FOR LASER RESPONSE
%save LaserResponse Spikes mergingPackage bads -append

%% looking at individual conditions and clusters for good clusters
%
plotRaster = true; % No raster plot in the figures!!
close all
for i=1:numel(Spikes)
    if ~ismember(i,bads)

        for I=4 %pick out conditions to look at
            ppms=fs/1000;
            spikes=(Spikes{i})*1000*ppms; %spikes back in samples
            name=Names{i};
            timeBefore=5000*ppms;timeAfter=8200*ppms;
            plotit=1;binsize=100*ppms;
            triggers=Conditions{I}.Triggers;
            [sp h bins trig_mech trig_light f]=...
                triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit,plotRaster);
            title(['Cluser ID: ', Names{i}, ', #: ',num2str(i)])
        end
        
    end
end


%% looking at possible artifacts


plotRaster = true; % raster plot in the figures!!
close all

for i=[1 2 3 4 5 6] % insert cluster numbers here
    for I=1 %pick out conditions to look at
        ppms=fs/1000;
        spikes=(Spikes{i})*1000*ppms; %spikes back in samples
        name=Names{i};
        timeBefore=5000*ppms;timeAfter=8200*ppms;
        plotit=1;binsize=10*ppms;
        triggers=Conditions{I}.Triggers;
        [sp h bins trig_mech trig_light f]=...
            triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit,plotRaster);
        title(['Cluster: ',num2str(i)])
    end
    %AM: repeated block to be able to look at all conditions simultaneuosly
    for I=2 %pick out conditions to look at
        ppms=fs/1000;
        spikes=(Spikes{i})*1000*ppms; %spikes back in samples
        name=Names{i};
        timeBefore=5000*ppms;timeAfter=8200*ppms;
        plotit=1;binsize=10*ppms;
        triggers=Conditions{I}.Triggers;
        [sp h bins trig_mech trig_light f]=...
            triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit,plotRaster);
        title(['Cluster: ',num2str(i)])
    end
    for I=3 %pick out conditions to look at
        ppms=fs/1000;
        spikes=(Spikes{i})*1000*ppms; %spikes back in samples
        name=Names{i};
        timeBefore=5000*ppms;timeAfter=8200*ppms;
        plotit=1;binsize=10*ppms;
        triggers=Conditions{I}.Triggers;
        [sp h bins trig_mech trig_light f]=...
            triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit,plotRaster);
        title(['Cluster: ',num2str(i)])
    end
    for I=4 %pick out conditions to look at
        ppms=fs/1000;
        spikes=(Spikes{i})*1000*ppms; %spikes back in samples
        name=Names{i};
        timeBefore=5000*ppms;timeAfter=8200*ppms;
        plotit=1;binsize=10*ppms;
        triggers=Conditions{I}.Triggers;
        [sp h bins trig_mech trig_light f]=...
            triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit,plotRaster);
        title(['Cluster: ',num2str(i)])
    end
end
%% Spikes backup.
% Elimination of the 'bads' in the 'Spikes' variable.
% But also making a backup before deleting. 
% WARNING! If you run this section twice without restoring the backup, you
% will lose good units!!
Spikes_BACKUP = Spikes;
Spikes(bads) = [];
%% Restore Spikes
% Use this section to restore the 'bads' labelled units into the 'Spikes'
% variable
Spikes = Spikes_BACKUP;
%% organize real figure


% population histogram for later plotting, by condition
ppms=fs/1e3;
spikes=cell2mat(Spikes)*1000*ppms;
plotit=1;
binsize=50*ppms;
timeBefore=2000*ppms;
timeAfter=8500*ppms;
H=[];
conds={};
count=0;
Trig_mech={};Trig_light={};
for I=1:numel(Conditions)   
    count=count+1;
    %spikes back in samples
    name=Names{I};
    triggers=Conditions{I}.Triggers;
    conds{count}=Conditions{I}.name;
    [sp, h, bins, Trig_mech{count}, Trig_light{count}, f]=...
        triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore,...
        timeAfter,Conditions{I}.name,mech,light,plotit);
    title(Conditions{I}.name)
    %convert to rate in Hz
    h=h*(1000/binsize*ppms);
    H(count,:)=h;
end
t=[-timeBefore:timeAfter]/ppms;
%% MORE BADS?
bads = [bads, 16];
bads = unique(bads);

%% get individual responses by condition

Sp={};
% close all
count=0;
plotit=1;

for I=1:numel(Conditions)  %by condition
    count=count+1;
    triggers=Conditions{I}.Triggers;
    SPIKES={};
    for i=1:numel(Spikes)   %for each neuron
        if consIdxs(i)
            spikes=(Spikes{i})*1000*ppms; %spikes back in samples
            name=Names{i};    
            %use same binning etc as pop hist
            [SPIKES{i} h bins trig_mech trig_light f]=triggeredAnalysisMUA(spikes,ppms,triggers,binsize,timeBefore, timeAfter,Conditions{I}.name,mech,light,plotit);
            title(['Cluster: ',num2str(i)])
        end
    end
    Sp{count}=SPIKES;
end



%

%% get data appropriate for plotting rasters


SPIKESs={};YSs={}; %all data here ,cond x num neuron

for ii=1:numel(Sp) %pick one condition
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

%% get color for each neuron, just for plotting
colormap jet;
cmap=colormap;
n=floor(size(cmap,1)/numel(SPIKESs{1}));
colors=cmap(1:n:end,:);
colors = jet(numel(SPIKESs{1}));

%% plot it all
figure('Color',[1,1,1])
% yUpLimit = max(cell2mat(cellfun(@(x) (cellfun(@max,x(end),'UniformOutput',false)),YSs)));
Ncl = numel(SPIKES);
Ncon = numel(Conditions);
for ii=1:Ncon
    auxAx = subplot(6,Ncon,[ii, ii+Ncon, ii+(2*Ncon)]);
    Nt = numel(Conditions{ii}.Triggers);
    for j=1:numel(SPIKESs{ii})        
        xs=SPIKESs{ii}{j}/ppms;
        ys=YSs{ii}{j};
        plot(xs,ys,'.','color',colors(j,:),'markersize',10)
        hold on
    end
    lstWV = find(~cellfun(@isempty,YSs{ii}),1,'last');
    yUpLimit = max(YSs{ii}{lstWV});
    yticks = ((1:Ncl) - 0.5) * Nt;
    set(auxAx,'YTick',yticks,'YTickLabel',1:Ncl);ylabel(...
        sprintf('Clusters_{(t = %d)}',Nt))
    
    % ylabel 'trials/neuron'
    box off
    xlim([min(bins) max(bins)])
    ylim([0,yUpLimit])
    
end

yLimit = max(H(:)) * 1.05;
%load(fname,'chan21')
%[~,cStack] =...
%     getStacks(false(1,length(chan21)),...
%     Conditions{1}.Triggers,'on',[-min(bins),max(bins)]*1e-3,...
%     fs,fs,[],chan21');
% meanMech = mean(squeeze(cStack),2);
for ii=1:Ncon
    
    auxAx = subplot(6,Ncon,[ii+(4*Ncon), ii+(5*Ncon)]);
    bar(bins,H(ii,:),'k')
    if ~strcmpi(Conditions{ii}.name,'lasercontrol')
        
    end
    xlim([min(bins) max(bins)])
    xlabel ms
    % ylabel 'pooled spike rate'
    title(Conditions{ii}.name)
    box off
    ylim([0 yLimit])
end

%
t = 0:1/fs:(length(Trig_mech{1}(1,:))-1)/fs;
t = 1000*(t - timeBefore/fs);
subplot(6,Ncon,4*Ncon - 3)
plot(t,Trig_mech{1}(1,:),'Color',[255, 128, 0]/255,'linewidth',2)
hold on;plot(t,...
    scaleOnLine(meanMech(1:length(t)),0,1),'Color',[255, 51, 0]/255,'linewidth',2)
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])
ylabel Stimulus

subplot(6,Ncon,4*Ncon - 2 )
plot(t,Trig_mech{2}(1,:),'Color',[255, 128, 0]/255,'linewidth',2);
hold on;plot(t,...
    scaleOnLine(meanMech(1:length(t)),0,1),'Color',[255, 51, 0]/255,'linewidth',2)
hold on
plot(t,Trig_light{2}(1,:),'Color', [0, 64, 255]/255','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])

subplot(6,Ncon,4*Ncon - 1)
plot(t,Trig_mech{3}(1,:),'Color',[255, 128, 0]/255,'linewidth',2);
hold on;plot(t,...
    scaleOnLine(meanMech(1:length(t)),0,1),'Color',[255, 51, 0]/255,'linewidth',2)
hold on
plot(t,Trig_light{3}(1,:),'Color', [0, 64, 255]/255','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])

subplot(6,Ncon,Ncon*4)
plot(t,Trig_light{4}(1,:),'Color', [0, 64, 255]/255','linewidth',1);
set(gca,'Visible','off')
box off
xlim([min(bins) max(bins)])
ylim([0 1.5])

%% plot some raw data
figure
load('M47_C2_L6Mech+L61mWanalysis','filteredResponse')

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




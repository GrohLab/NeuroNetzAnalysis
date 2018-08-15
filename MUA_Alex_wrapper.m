%% import data from spike2 into mat
% fname='1300_1600_3320_whisker.smr'

fname='M137_C4_Mech+L6 05mW'; % EIC Sailaja's analysis
% fdir='C:\Users\neuro\Documents\MATLAB\16 channel\';
fdir = 'F:\Experiments_2018\19_4_2018';
%X=importSMR([fname,'.mat'],fdir,0);

%%
% 
% load ('whisker_lpf.mat');
% whsk=whisking_lpf;
% whsk=whsk/max(whsk)+1;
%% make input for clustering
clc
cd(fdir)
load([fname, '.mat'])

num_trials   = 1;
% Fs           = 20000; The sampling frequency is taken from the file.
Fs           = head1.SamplingFrequency;
num_channels =  16;
trial_dur    = 100;

% chanOrder=[8 9 7 10 4 13 5 12 2 15 1 16 6 11 3 14]; % poly design
chanOrder=[6 9 7 10 4 13 5 12 2 15 1 16 6 11 3 14]; % linear design
%chanOrder=fliplr(chanOrder);

Data={};
for i=1:num_trials
    for j = 1:num_channels
        eval(['ns(j)=numel(chan' num2str(j) ');']);j
    end
    
    
    ns=min(ns);
    data=zeros(ns,num_channels);
    for j = 1:num_channels
        eval(['temp=double(chan' num2str(j) ');']);
        data(:,j) = temp(1:ns); j
    end
    
    data=data(:,chanOrder);%reorder
    Data{i}=data;
    clear data;
    clear chan*
    clear head*
    clear temp
    clear ns
end

% save Data_1300_1600_3415_w_puff_marker_filtered.mat Data '-v7.3'
save([fname,'data.mat'], 'Data', '-v7.3')
disp('imported data into "Data"')


%% run the detection and clustering functions

% 1:4 3:6 7:10 12:15
% 1:4 3:6 7:10 5:8 9:12
channelPacks = {1:4 4:7 7:10 10:13 13:16};
chanPackIdx = 5;
data={};
ch = channelPacks{chanPackIdx}; %  channels for sorting
for i=1:numel(ch)
    data{1,1}(:,i)=Data{1,1}(:,ch(i));
end


% thresh= 3.9; % number of standard deviations above background noise
thresh= 5; % number of standard deviations above background noise
spikes = ss_default_params(Fs, 'thresh', thresh);
spikes = ss_detect(data,spikes);
spikes = ss_align(spikes);
spikes = ss_kmeans(spikes);
spikes = ss_energy(spikes);
spikes = ss_aggregate(spikes);

% main tool
splitmerge_tool(spikes)

%%
% stand alone outlier tool
outlier_tool(spikes)
%% Save
save([fname,'channel',num2str(chanPackIdx),'.mat'],'spikes')
%% collect clusters into cell array sortedData (exlude garbage clusters)

% fdir='C:\Users\alex\Desktop\16chanelanalysis\120213\01_whiskeronly_sorted\';
% cd(fdir);
str='Pack';exclude={};
chanData=matchfiles(filePath,str,exclude); % helper function collects folders
sortedData={};

for j=1:numel(chanData)
    load(chanData{j,1});
    ktemp=size(sortedData)+1; k=ktemp(1);
    name=chanData{j}; name=name(end-8:end-4);
    for i=1:length(spikes.labels)
        indices=find(spikes.assigns==spikes.labels(i,1));
        if ~isempty(indices) && spikes.labels(i,2)~=4;
           
            sortedData{k,1}=(['ch' name '_cl_' num2str(spikes.labels(i,1))]); %channel and cluster
            sortedData{k,2}=spikes.spiketimes(indices); %spikes times
            sortedData{k,3}=spikes.labels(i,2); %cluster labels: 2=good, 4=schrott
            k=k+1;
        end
    end
end


% include= [15: 20]; % [3 5 17]
%% sortedData=sortedData(include,:) % which units ?
% Avoid overwritting the sorted spikes for two data files in the same
% folder
save([fname,'_all_channels.mat'],'sortedData')
% run till here and then move to RM script

%% raster plot and coincident spike times

d=0.002; %window for coincidence detection (d before spike)
% d=1/10000;
includeIdenticalSp=0; % include identical spike times in column 4 (i.e. double sorted spike times)
time_a=405; %only for raster plot
time_b=430; %only for raster plot
overlap=5; % thresh for indicating overlapping spiketimes in %
samplRate=20000;

tSpike=time_a:1/samplRate:time_b-1/samplRate;

% figure(1)
% subplot(6,1,1);plot(tSpike, whsk(time_a*samplRate+1:time_b*samplRate), 'k')
% subplot(6,1,2:6);
[sortedData, nCSP]=coincSpikes(sortedData, d, time_a, time_b, overlap, includeIdenticalSp, tSpike, samplRate, whsk);



figure, bar3(nCSP)
figure, pcolor(nCSP), colormap(gray(10)), shading flat, colorbar

% save sortedData_allclusters.mat sortedData '-v7.3'

% content of sortedData

% 1channelsAndClustername 2spiketimes 3label 4coincidentspikes 5clusternameConnectedTo 6percetageOfCoincidentSpikes 7convolvedSpiketrain 

%% make histograms
allCoinc=[];
allSpikes=[];
bin=0.005;
for i= 1:length(sortedData)
    allCoinc=[allCoinc sortedData{i,4}];
    allSpikes=[allSpikes sortedData{i,2}];
end

x=[time_a:bin:time_b];

hCoinc=histc(allCoinc, x); % hCoinc=smooth(hCoinc/max(hCoinc));
hAll=histc(allSpikes, x); % hAll=smooth(hAll/max(hAll));

figure
bar(x, hCoinc, 'r')
 hold on
 plot(tSpike, whsk(time_a*samplRate+1:time_b*samplRate), 'k')
 hold off
axis tight
figure
bar(x, hAll, 'k') 
axis tight
%%
% main tool
splitmerge_tool(spikes)
%%
% Note: In the code below, "clus", "clus1", "clus2", and "clus_list" are dummy
% variables.  The user should fill in these vlaues with cluster IDs found
% in the SPIKES object after running the algorithm above.
%

% plots for single clusters
plot_waveforms( spikes, clus );
plot_stability( spikes, clus);
plot_residuals( spikes,clus);
plot_isi( spikes, clus );
plot_detection_criterion( spikes, clus );

% comparison plots
plot_fld( spikes,clus1,clus2);
plot_xcorr( spikes, clus1, clus2 );

% whole data plots
plot_features(spikes );
plot_aggtree(spikes);
show_clusters(spikes, [clus_list]);
compare_clusters(spikes, [clus_list]);

% outlier manipulation (see M-files for description on how to use)
spikes = remove_outliers( spikes, which );
spikes = reintegrate_outliers( spikes, indices, mini );

% quality metric functions
%
% Note: There are versions of these functions in the quality_measures
% directory that have an interface that does not depend on the SPIKES
% structure.  These are for use by people who only want to use the quality
% metrics but do not want to use the rest of the package for sorting.
% These functions have the same names as below but without the "ss_" prefix.
%
FN1 = ss_censored( spikes, clus1 )  % 0.0015 juxta, 0 for cluster containing different spike parts from same neuron
FP1 = ss_rpv_contamination( spikes, clus1  ) % 0 juxta, 1 for cluster containing different spike parts from same neuron
FN2 = ss_undetected(spikes,clus1) % 1.1768e-039 juxta, 0.2563 for cluster containing different spike parts from same neuron
confusion_matrix = ss_gaussian_overlap( spikes, clus1, clus2 )




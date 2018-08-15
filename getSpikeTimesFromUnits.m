function [sortedData] = getSpikeTimesFromUnits(rootNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% str='_Pack';
% chanData=matchfiles(fdir,str,exclude); % helper function collects folders
% sortedData={};
sortedFile = fullfile([rootNames,'_all_channels.mat']);
if exist(sortedFile,'file')
    disp('The file exist')
    load(sortedFile,'sortedData')
    anw = input('Would you like to overwrite it? (y/n) ','s');
    if anw~='y'
        return
    end    
end
chanData = dir([rootNames,'_Pack*.mat']);
auxData = cell(1,numel(chanData));
catFields = @(x) fullfile(x.folder,x.name);
for cf = 1:numel(chanData)
    auxData(cf) = {catFields(chanData(cf))};
end
chanData = auxData;
sortedData = {};
for j=1:numel(chanData)
    load(chanData{j,1});
    ktemp=size(sortedData)+1; k=ktemp(1);
    name=chanData{j}; name=name(end-8:end-4);
    for i=1:length(spikes.labels)
        indices=find(spikes.assigns==spikes.labels(i,1));
        if ~isempty(indices) && spikes.labels(i,2)~=4
            sortedData{k,1}=(['ch' name '_cl_' num2str(spikes.labels(i,1))]); %channel and cluster
            sortedData{k,2}=spikes.spiketimes(indices); %spikes times
            sortedData{k,3}=spikes.labels(i,2); %cluster labels: 2=good, 4=schrott
            k=k+1;
        end
    end
end

save(sortedFile,'sortedData')
end
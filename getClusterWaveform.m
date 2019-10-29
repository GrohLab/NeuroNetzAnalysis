function waveform = getClusterWaveform(clusterID, dataDir)
%GETCLUSTERWAVEFORM reads the raw binary file and compiles the voltage
%traces for the given cluster (the cluster number should be the one
%assigned by Kilosort/Phy). The output is a cell (or a structure)
%containing the mean waveform and its standard deviation.
%   waveform = getClusterWaveform(clusterID)
%       INPUTS
%           - clusterID - character array, double or cell array containing
%           the clusters from which the waveforms are required.
%       OUTPUT
%           - waveform - cell array or double vector containing the mean
%           wavefrom for the given cluster(s).
% Emilio Isaias-Camacho @GrohLab 2019

%% Input validation
waveform = [];
checkNature = @(x) [iscell(x), ischar(x), isnumeric(x)];
if ~any(checkNature(clusterID))
    fprintf(1,'Unsupported input!\n')
    fprintf(1,'Be sure you input either the cluster ID as a ')
    fprintf(1,'character vector, a number or\nas a cell array')
    return
end
if ~exist(dataDir, 'dir')
    fprintf(1,'Not possible to retrieve waveforms without the data!\n')
    fprintf(1,'Please provide the data directory\n')
    return
end
%% Getting ready for the file reading
% Reading the cluster summary
clTable = readClusterInfo(fullfile(dataDir, 'cluster_info.tsv'));
% Converting the ID(s) to cell arrays according to their nature
switch bi2de(checkNature(clusterID),'left-msb')
    case 1
        fprintf(1,'Numeric ID detected\n')
        clusterID = {num2str(clusterID)};
    case 2
        fprintf(1,'Character ID detected\n')
        clusterID = {clusterID};
    case 4
        fprintf(1,'Cell ID detected\n')
        charFlag = all(cellfun(@ischar, clusterID));
        if ~charFlag
            fprintf(1,'The cluster ID(s) should be only character within')
            fprintf(1,' a cell array\n')
            fprintf(1,'Please provide the ID(s) as required and try again\n')
            return
        end
end
% Determining hosting channels
try
    ch2read = clTable{clusterID, 'channel'};
catch
    fprintf(1,'Some of the given clusters do not exist in this experiment\n')
    fprintf(1,'Clusters not found:\n')
    missClustFlag = false(numel(clusterID),1);
    clIdx = false(numel(clTable.Properties.RowNames), numel(clusterID));
    for ccl = 1:numel(clusterID)
        clIdx(:,ccl) = strcmp(clTable.id, clusterID(ccl));
        if ~any(clIdx(:,ccl))
            missClustFlag(ccl) = true;
            fprintf(1,'%s\n', clusterID{ccl})
        end
    end
    if sum(missClustFlag) < numel(clusterID)
        contAns = questdlg('Continue without these clusters?', 'Continue?',...
            'Yes','No','Yes');
        if strcmp(contAns, 'No')
            fprintf(1,'Aborting...\n')
            return
        end
    else
        fprintf(1,'No valid cluster ID provided!\n')
        return
    end
    clusterID(missClustFlag) = [];
    clIdx(:,missClustFlag) = [];
    ch2read = clTable{clusterID, 'channel'};
end

% Determinig the spike times for the given clusters
spikeFile = dir(fullfile(dataDir,'*_all_channels.mat'));
if ~isempty(spikeFile)
    load(fullfile(dataDir, spikeFile.name), 'sortedData', 'fs')
    if ~exist('fs','var')
        fsFile = dir(fullfile(dataDir,'*_sampling_frequency.mat'));
        load(fullfile(dataDir, fsFile.name), 'fs')
    end
end
%% Reading the binary file
binFile = dir(fullfile(dataDir, '*.bin'));
if isempty(binFile)
    fprintf(1, 'Without a binary file it is impossible to get the waveforms')
    fprintf(1,'\n');
    return
end
spkSubs = cellfun(@(x) round(x.*fs),sortedData(:,2),...
    'UniformOutput',false);
end


function waveTable = getClusterWaveform(clusterID, dataDir)
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
waveTable = table();
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
clusterID = unique(clusterID);
%% Getting ready for the file reading
% Reading the cluster summary
clTable = readClusterInfo(fullfile(dataDir, 'cluster_info.tsv'));
% Reading the channel map
chanMap = readNPY(fullfile(dataDir, 'channel_map.npy'));
Nch = numel(chanMap)-1;
% Verifying if the given clusters exist in this experiment
clIdx = false(size(clTable, 1), numel(clusterID));
for ccl = 1:numel(clusterID)
    clIdx(:,ccl) = strcmp(clTable.id, clusterID(ccl));
end
missClustFlag = ~any(clIdx,1);
if ~all(~missClustFlag)    
    fprintf(1,'Some of the given clusters do not exist in this experiment\n')
    fprintf(1,'Clusters not found:\n')
    fprintf(1,'%s\n', clusterID{missClustFlag})
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
end    
clusterID(missClustFlag) = [];
clIdx(:,missClustFlag) = [];
clIdx = any(clIdx,2);
% Determining hosting channels
ch2read = chanMap(clTable{clusterID, 'channel'} + 1);


% Verifying if the waveform(s) for the given cluster(s) was/were computed
% already
%{
waveFile = dir(fullfile(dataDir,'_waveforms.mat'));
if exist(waveFile, 'file')
    load(fullfile(dataDir, waveFile.name),'waveTable')
    N_exCl = size(waveTable, 1);
    exIdx = false(N_exCl, numel(clusterID));
    for ccl = 1:N_exCl
        exIdx(:,ccl) = strcmp(waveTable.id, clusterID(ccl));
    end
end
%}

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
% Taking ?1.25 ms around the spike.
spikeWaveTime = 2*round(1.25e-3 * fs) + 1;
spikeSamples = (spikeWaveTime - 1)/2;
binFile = dir(fullfile(dataDir, '*.bin'));
if isempty(binFile)
    fprintf(1, 'Without a binary file it is impossible to get the waveforms')
    fprintf(1,'\n');
    return
end
spkSubs = cellfun(@(x) round(x.*fs),sortedData(clIdx,2),...
    'UniformOutput',false);
waveform = zeros(numel(clusterID), spikeWaveTime);
% [ch2read, readOrder, repeatChs] = unique(ch2read);
fID = fopen(fullfile(dataDir, binFile.name), 'r');
cchan = 1;
% Main loop
while ~feof(fID) && cchan <= numel(clusterID)
    fprintf(1,'Reading channel %d ',ch2read(cchan))
    % Jumping to the channel 
    fseek(fID, 2*(ch2read(cchan)), 'cof');
    spkDists = [spkSubs{cchan}(1);diff(spkSubs{cchan})];
    fprintf(1,'looking for cluster %s...', clusterID{cchan})
    for cspk = 1:numel(spkSubs{cchan})
        % Jumping to 1 ms before the time when the spike occured
        fseek(fID, 2*((Nch+1)*(spkDists(cspk) - spikeSamples)), 'cof');
        % Reading the waveform
        rawWave = fread(fID, [1, spikeWaveTime], 'int16=>single', 2*Nch);
        % Averaging the waveform
        waveform(cchan,:) = waveform(cchan, :) + rawWave /...
            cspk;
        % Jumping back to the exact time of the spike
        fseek(fID, -2*((Nch+1)*(spikeSamples+1)), 'cof');
    end
    fprintf(1,' done!\n')
    cchan = cchan + 1;
    frewind(fID);
end
fclose(fID);

%% Arranging the output
if isrow(clusterID)
    clusterID = clusterID';
end
waveTable = table(waveform, ch2read,...
    'RowNames', clusterID, 'VariableNames', {'Waveform', 'Channel'});
waveTable.Properties.DimensionNames{1} = 'id';
end


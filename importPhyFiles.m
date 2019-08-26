function importPhyFiles(pathToData, outputName, removeNoise, ChAndAmpFlag)
%IMPORTPHYFILES reads the files created by both Kilosort and Phy or Phy2
%and, depending on the flag configurations, it saves a .mat file. The
%contents of this file are constructed with the output files (.tsv, .npy).
%           importPhyFiles( pathToData, outputName, removeNoise,
%               ChAndAmpFlag)
%       INPUTS: no marking means required, [] mean optional but needed for
%                                             later input arguments.
%           pathToData - Path string to where the .tsv and .npy files are
%                           located
%           [outputName] - Name for the .mat file to be created.
%                           WARNING! if the file already exists, the
%                           function will overwrite it
%           [removeNoise] - Default false; flag to indicate the inclusion
%                           of the clusters labelled as noise.
%           [ChAndAmpFlag] - Default false; flag to indicate the inclusion
%                            of channel and amplitude of each cluster.
%       OUTPUTS:
%           No variable imported to the workspace but a .mat files written
%           in the given folder; either in the pathToData or in the
%           outputName.
%  Alexandra Merkel & Emilio Isaias-Camacho @ GrohLab 2019
% Maximum number of groups. 3 so far: good, MUA, and noise

MX_CLSS = 3;
% Take Kilosort Data and create *all_channels.map for further analysis
addpath(genpath('C:\Users\NeuroNetz\Documents\npy-matlab')) % path to npy-matlab
addpath(genpath('C:\Users\NeuroNetz\Documents\GitHub\NeuroNetzAnalysis')) % path to readTSV

% Where is your KiloSort Data saved
%pathToData = 'C:\Users\NeuroNetz\Documents\Data\Alex\M73_attempt3';

% What is your outputfile called
if ~exist('outputName','var') || isempty(outputName)
    binFile = dir(fullfile(pathToData,'*.bin'));
    binName = binFile(1).name;[~,bfName,fExt] = fileparts(binName);
    outputName = strrep(bfName,fExt,'');
end
% Input sampling frequency (not written in the phy files from KiloSort)
try
    load(fullfile( pathToData, [outputName, '_sampling_frequency.mat']),...
        'fs')
catch LE
    disp( LE.message)
    fsString = inputdlg('Please provide a sampling frequency:');
    try
        fs = str2double(fsString{1});
        if isnan(fs)
            fprintf('I''m sorry... you should put in ONLY NUMBERS :)\nStart again\n')
            fprintf('No output file written\n')
            return
        end
    catch 
        fprintf(1,'Cancel button pressed. ')
        fprintf(1,'Transforming the files might help you get the sampling ')
        fprintf(1,'frequency.\n');
        return
    end
    save(fullfile( pathToData, [outputName, '_sampling_frequency.mat']),...
        'fs')
end

% Do you want to remove Noise Data?
if ~exist('removeNoise','var')
    removeNoise = false;
end

% Do you want to include channel and amplitude information?
if ~exist('ChAndAmpFlag','var')
    ChAndAmpFlag = false;
end

% Reading the KiloSort datafiles (phy updates these files)
spkTms = readNPY(fullfile(pathToData, 'spike_times.npy'));
clID = readNPY(fullfile(pathToData, 'spike_clusters.npy'));
clGr = readTSV(fullfile(pathToData,'cluster_group.tsv'));
if ChAndAmpFlag
    clInfo = getClusterInfo(fullfile(pathToData,'cluster_info.tsv'));
end

nIdx = cellfun(@strcmp,clGr(:,2),repmat("noise",size(clGr,1),1));
gIdx = cellfun(@strcmp,clGr(:,2),repmat("good",size(clGr,1),1));
mIdx = cellfun(@strcmp,clGr(:,2),repmat("mua",size(clGr,1),1));
grIdx = [gIdx, mIdx, nIdx];

allClusters = cellfun(@num2str,clGr(:,1),'UniformOutput',false);

spkCell = cell(size(clGr,1),1);

for ccln = 1:numel(allClusters)
    spkCell(ccln) = {double(spkTms(clID == clGr{ccln,1}))/fs};
end

row = find(grIdx');
[r,col2] = ind2sub(size(grIdx),row);
[~,rearrange] = sort(r);
clGroup = col2(rearrange);
[row2, ~] = find(grIdx);
cuts = sum(grIdx,1);
limits = [0,cumsum(cuts)];
for cgroup = 1:MX_CLSS
    clGroup(row2(limits(cgroup)+1:limits(cgroup+1))) = cgroup;
end
sortedData = cat(2,allClusters,spkCell,num2cell(clGroup));
if ChAndAmpFlag
    importedFlag = ismember(clInfo.id,sortedData(:,1));
    CHAN_AMPS = [clInfo.channel(importedFlag), clInfo.Amplitude(importedFlag)];
else
    CHAN_AMPS = [];
end
% Removes the noise from the matrix and saves an alternative output file
if removeNoise
    index = cellfun(@(x) x==3,sortedData(:,3));
    sortedData(index,:) = []; %#ok<NASGU>
    if ChAndAmpFlag
        CHAN_AMPS(index,:) = []; %#ok<NASGU>
    end
    filename = [outputName, '_all_channels.mat'];
    fname = fullfile(pathToData, filename);
    save(fname, 'sortedData', 'fs','CHAN_AMPS');
else
    filename = [outputName, '_all_channels.mat'];
    fname = fullfile(pathToData, filename);
    save(fname, 'sortedData', 'fs','CHAN_AMPS');
end
end
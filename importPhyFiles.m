function [ sortedData ] = ...
    importPhyFiles(pathToData, outputName, removeNoise)
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
%       OUTPUTS:
%           No variable imported to the workspace but a .mat files written
%           in the given folder; either in the pathToData or in the
%           outputName.
%  Alexandra Merkel & Emilio Isaias-Camacho @ GrohLab 2019
% Maximum number of groups. 4 so far: 1 - good, 2 - MUA, 3 - noise, and 4 -
% <undefined>

grps = ["good", "mua", "noise", ""];
% Take Kilosort Data and create *all_channels.map for further analysis

fnOpts = {'UniformOutput', false};

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

% Reading the KiloSort datafiles (phy updates these files)
spkTms = readNPY(fullfile(pathToData, 'spike_times.npy'));
clID = readNPY(fullfile(pathToData, 'spike_clusters.npy'));
clGr = readTSV(fullfile(pathToData,'cluster_group.tsv'));

% Creating a logical index for the user labels
[~, clGroup] = cellfun(@(x) ismember(x, grps), clGr(:,2));
allClusters = cellfun(@num2str,clGr(:,1), fnOpts{:});
spkCell = cellfun(@(x) double(spkTms(clID == x))/fs, clGr(:,1), fnOpts{:});
sortedData = cat(2,allClusters,spkCell,num2cell(clGroup));
% Removes the noise from the matrix and saves an alternative output file
if removeNoise
    sortedData(nIdx,:) = [];
end
filename = [outputName, '_all_channels.mat'];
fname = fullfile(pathToData, filename);
if exist(fname, 'file')
    saveAns = questdlg(sprintf('%s exists already! Overwrite?',filename),...
        'Overwrite?', 'Yes', 'No', 'Yes');
    if strcmp(saveAns, 'No')
        fprintf(1, 'No file written!\n')
        return
    end
end
save(fname, 'sortedData', 'fs');
end
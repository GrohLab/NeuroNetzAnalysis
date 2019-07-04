function importPhyFiles(pathToData, outputName, removeNoise)

% Take Kilosort Data and create *all_channels.map for further analysis
addpath(genpath('C:\Users\NeuroNetz\Documents\npy-matlab')) % path to npy-matlab
addpath(genpath('C:\Users\NeuroNetz\Documents\GitHub\NeuroNetzAnalysis')) % path to readTSV

% Where is your KiloSort Data saved
%pathToData = 'C:\Users\NeuroNetz\Documents\Data\Alex\M73_attempt3';

% What is your outputfile called
if ~exist('outputName','var')
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
    
end

% Do you want to remove Noise Data?
if ~exist('removeNoise','var')
    removeNoise = false;
end

% Reading the KiloSort datafiles (phy updates these files)
spkTms = readNPY(fullfile(pathToData, 'spike_times.npy'));
clID = readNPY(fullfile(pathToData, 'spike_clusters.npy'));
clGr = readTSV(fullfile(pathToData,'cluster_group.tsv'));

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
limits = cumsum(cuts);
clGroup(row2(1:limits(1))) = 1;
clGroup(row2(limits(1)+1:limits(2))) = 2;
clGroup(row2(limits(2)+1:limits(3))) = 3;
sortedData = cat(2,allClusters,spkCell,num2cell(clGroup));

% Removes the noise from the matrix and saves an alternative output file
if removeNoise
    index = sortedData{:,3} == 3;
    sortedData(index,:) = []; %#ok<NASGU>
    filename = [outputName, '_all_channels.mat'];
    fname = fullfile(pathToData, filename);
    save(fname, 'sortedData', 'fs');
else
    filename = [outputName, '_all_channels.mat'];
    fname = fullfile(pathToData, filename);
    save(fname, 'sortedData', 'fs');
end
end
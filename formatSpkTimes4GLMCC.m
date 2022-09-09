function [iOk] = formatSpkTimes4GLMCC(spkTms, dataDir, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser();
checkST = @(x) iscell(x) & all(cellfun(@(y) isnumeric(y), x)) & ...
    isvector(x);
checkDD = @(x) (ischar(x) | isstring(x)) & exist(x, 'dir');
p.addRequired('spkTms', checkST)
p.addRequired('dataDir', checkDD)

p.parse(spkTms, dataDir, varargin{:})

spkTms = p.Results.spkTms;
dataDir = p.Results.dataDir;
glmccDir = fullfile(dataDir, "GLMCC");

if ~exist(glmccDir, 'dir')
    if ~mkdir(glmccDir)
        fprintf(1, 'Problem while creating directory\n')
        fprintf(1, 'Keeping the given directory as output directory!\n')
        % User input; continue or abort; give another directory;
        glmccDir = dataDir;
    end
else
    cellFiles = dir(glmccDir);
    cellFiles([cellFiles.isdir]) = [];
    % User input: overwirte or abort; choose cells;
end

if ~exist('cellFiles','var') || isempty(cellFiles) && ...
        numel(cellFiles) < size(spkTms,1)
    for ccl = 1:size(spkTms,1)
        fID = fopen(fullfile(glmccDir, sprintf('cell%d.txt', ccl-1)), 'w');
        if fID >= 3
            fprintf(fID, "%f\n", spkTms{ccl});
        end
        fclose(fID);
    end
end

all_path = fullfile(glmccDir,'all.txt');
if ~exist(all_path, 'file')
    fID = fopen(all_path, "w");
    cellfun(@(x) writeSpikesFile(fID, x), spkTms);
    fclose(fID);
end
end

function writeSpikesFile(fID, spkTms)
    fprintf(fID, "%f\n", spkTms);
    fprintf(fID, ";\n");
end
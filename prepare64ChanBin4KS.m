function iOk = prepare64ChanBin4KS(dataDir, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Input parsing
p = inputParser;

checkDataDir = @(x) exist(x,'dir');

defOutBinName = 'MedianFiltered_int16';
checkOutBinName = @(x) all([~isempty(x), ischar(x) | isstring(x)]);

defFs = 3e4;
checkFs = @(x) all([isnumeric(x), x > 0, numel(x) == 1]);

defChan = 64;
checkChan = @(x) all([isPositiveIntegerValuedNumeric(x), numel(x) == 1]);

p.addRequired('dataDir', checkDataDir)
p.addParameter('BinFileName', defOutBinName, checkOutBinName);
p.addParameter('fs', defFs, checkFs);
p.addParameter('ChanNumber', defChan, checkChan);

parse(p, dataDir, varargin{:});

dataDir = p.Results.dataDir;
outBinName = p.Results.BinFileName;
fs = p.Results.fs;
Nch = p.Results.ChanNumber;

%% Validation of inputs
iOk = false; getName = @(x) string(x.name);
% Search for Recording*.bin files in the given folder
recFiles = dir(fullfile(dataDir, 'Recording*.bin'));
if isempty(recFiles)
    fprintf(1, 'No recording file found!\n')
    fprintf(1, 'Make sure you provide a directory with a recording file');
    fprintf(1, ' in it!\n');
    return
end

% Validation for existing output file
outBinName = string(outBinName);
if ~contains(outBinName, ".bin")
    outFullName = fullfile(dataDir, outBinName + ".bin");
else
    outFullName = fullfile(dataDir, outBinName);
end
if exist(outFullName,'file')
    ovwtAns = questdlg(...
        sprintf('Warning! File %s exists. Overwrite?',outFullName),...
        'Overwrite?','Yes','No','No');
    if strcmpi(ovwtAns,'no')
        fprintf('No file written.\n')
        return
    end
end
%% Processing files
% How many files are in the folder
Nrf = numel(recFiles);
% If there are more than 1 files in the folder, ask the used which are
% going to be merged.
if Nrf > 1
    selFile = listdlg("ListString", arrayfun(getName, recFiles));
    % Of course, if the user presses cancel, the function ends without
    % creating any file.
    if isempty(selFile)
        fprintf(1, "Canceled. No file created!\n")
        return
    end
    % Getting only those desired files.
    recFiles = recFiles(selFile); Nrf = numel(recFiles);
    fNames = arrayfun(getName, recFiles);
    if Nrf > 1
        % Merging in a given order
        fileOrder = (1:Nrf)';
        defInput = num2cell(num2str(fileOrder));
        answr = inputdlg(fNames,'File order',[1, 60], defInput);
        nFileOrder = str2double(answr);
        if ~isempty(answr) && sum(abs(fileOrder - nFileOrder)) ~= 0
            fprintf(1,'File order changed...\n')
            recFiles = recFiles(nFileOrder);
        end
    end
end
fNames = arrayfun(getName, recFiles);
try
    mem = memory;
catch
   fprintf(1, "Memory pack unavailable!\n")
end
% Overwrite the output file.
fPerm = 'w';
for cf = 1:Nrf
    fprintf(1, "File: %s\n", fNames(cf));
    % How many bytes is this file?
    fWeight = recFiles.bytes;
    % How many bytes are available in this computer. Occupying 85%
    mxBytes = 0.8 * mem.MaxPossibleArrayBytes;
    % Number of samples in the file with Nch channels. 2 bytes from uint16.
    % 4 from int32 and considering the samples that fit in memory. 
    Ns = fWeight/(2*Nch); Nms = mxBytes/(4*Nch);
    if Ns > Nms
        % If the file has more samples than the memory allows, then take
        % the samples that fit in memory.
        Ns = Nms;
    else
        % If the file fits entirely in the memory, then read it all at once
        Ns = Inf;
    end
    % Low-level programming keeps the current position in the opened data
    % file.
    dfID = openDataFile(fullfile(dataDir, fNames(cf)));
    chnk = 1;
    while ~feof(dfID)
        fprintf(1, "Reading batch %d...\n", chnk)
        chnk = chnk + 1;
        dataBuff = readDataFileAndMedianFilter(dfID, Ns, Nch);
        ofID = createOutputFile(outFullName, fPerm);
        fprintf(1, "Writing...\n")
        fwrite(ofID, dataBuff, "int16")
        fPerm = 'a'; fclose(ofID);
    end
    fclose(dfID);
end
iOk = true;
%% Writing useful files
ffoID = fopen(fullfile(dataDir,outBinName + "_fileOrder.txt"),'w');
fprintf(ffoID,"%s\n",arrayfun(getName, recFiles));
fprintf(ffoID, "%s.bin", outBinName); fclose(ffoID);
save(fullfile(dataDir, outBinName + "_sampling_frequency.mat"), 'fs')
end

function fID = openDataFile(fPath)
fID = fopen(fPath, 'r');
if fID < 3
    % Error by opening data file
    [~, fName] = fileparts(fPath);
    fprintf(1, 'Unable to open %s!\n', fName)
    fprintf(1, 'Verify that no other process is running with it and')
    fprintf(1, ' try again.\n');
    return
end
end

function [dataBuff, buffMed] = readDataFileAndMedianFilter(fID, Ne, Nch)
% As a result to removing the median of the channels, the data will take
% negative values. This is why I need to convert the data to signed 32-bit
% integers before substracting the median. Once the data is median  free,
% then it can be converted back to signed 16-bit integer.
dataBuff = fread(fID, [Nch, Ne], 'uint16=>int32');
buffMed = median(dataBuff);
dataBuff = int16(dataBuff - buffMed);
end

function fID = createOutputFile(fPath,fPerm)
fID = fopen(fPath, fPerm);
if fID < 3
    fprintf(1, 'Unable to create/open output file!\n');
    fprintf(1, 'Please verify the file and directory!\n');
    return
end % Data file identifier validation
end

function fID = createMedianFile(fPath) %#ok<DEFNU>
fID = fopen(fPath, 'w');
if fID < 3
    fprintf(1, 'Unable to create median file!\n');
    fprintf(1, 'Please verify the file and directory!\n');
    return
end % Median File identifier validation
end

%{
if numel(recFile) < 2
    % One recording file in directory
    dfID = openDataFile(fullfile(dataDir, recFile.name));
    memObj = memory;
    if 0.85*memObj.MaxPossibleArrayBytes > 2*recFile.bytes
        % It is possible to allocate the whole recording.
        [data, expMedian] = readDataFileAndMedianFilter(dfID, Inf);
        fclose(dfID);
        ofID = createOutputFile(fullfile(dataDir, [outBinName, '.bin']));
        fwrite(ofID, data, 'int16'); fclose(ofID);
        mfID = createMedianFile(fullfile(dataDir, 'ExperimentMedian.bin'));
        fwrite(mfID, expMedian, 'int16'); fclose(mfID); fs = 3e4; 
        save(fullfile(dataDir, [outBinName, '_sampling_frequency.mat']),...
            'fs')
        iOk = true;
    else
        % Read a big chunk of data, and write in disk. Use while loop.
        % TODO: Implement this functionality
        
    end % Memory validation
else
    % Deal with more than one recording file
    % TODO: Allow the user to select one recording file from the folder.
    fprintf(1, 'There''s more than 1 recording in the folder!/n')
    % Choosing files. If the user selects more than one, the function will
    % ask in which order to merge.
    recFileSubs = listdlg('ListString', arrayfun(@(x) x.name, recFile));
    Nrf = numel(recFileSubs);
    if Nrf > 1
    end
end % More than one Recording* file in the directory

%}
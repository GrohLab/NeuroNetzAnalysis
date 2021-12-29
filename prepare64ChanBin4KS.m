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
p.addOptional('BinFileName', defOutBinName, checkOutBinName);
p.addOptional('fs', defFs, checkFs);
p.addOptional('ChanNumber', defChan, checkChan);

p.parse(dataDir, varargin{:});

dataDir = p.Results.dataDir;
outBinName = p.Results.BinFileName;
fs = p.Results.fs;
Nch = p.Results.ChanNumber;

%% Validation of inputs
iOk = false;

recFile = dir(fullfile(dataDir, 'Recording*.bin'));
if isempty(recFile)
    fprintf(2, 'No recording file found!\n')
    fprintf(2, 'Make sure you provide a directory with a recording file');
    fprintf(2, ' in it!\n');
    return
end
% Validation for existing output file
outFullName = fullfile(dataDir, outBinName);
if exist(outFullName,'file')
    ovwtAns = questdlg(...
        sprintf('Warning! File %s exists. Overwrite?',outFullName),...
        'Overwrite?','Yes','No','No');
    if strcmpi(ovwtAns,'no')
        fprintf('No file written.\n')
        return
    end
end

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

function [dataBuff, buffMed] = readDataFileAndMedianFilter(fID, Ne)
dataBuff = fread(fID, Ne, 'uint16=>int32');
dataBuff = reshape(dataBuff, 64, []);
buffMed = int16(median(dataBuff));
dataBuff = int16(dataBuff - buffMed);
end

function fID = createOutputFile(fPath)
fID = fopen(fPath, 'w');
if fID < 3
    fprintf(1, 'Unable to create output file!\n');
    fprintf(1, 'Please verify the file and directory!\n');
    return
end % Data file identifier validation
end

function fID = createMedianFile(fPath)
fID = fopen(fPath, 'w');
if fID < 3
    fprintf(1, 'Unable to create median file!\n');
    fprintf(1, 'Please verify the file and directory!\n');
    return
end % Median File identifier validation
end
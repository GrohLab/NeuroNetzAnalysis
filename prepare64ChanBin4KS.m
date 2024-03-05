function iOk = prepare64ChanBin4KS(dataDir, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Input parsing
p = inputParser;

checkDataDir = @(x) exist(x,'dir');

defOutBinName = 'MedianFiltered_int16';
checkOutBinName = @(x) all([~isempty(x), ischar(x) | isstring(x)]);

defFs = 3e4;
checkFs = @(x) isnumeric(x) && (x > 0) && isscalar(x);

defChan = 64;
checkChan = @(x) isPositiveIntegerValuedNumeric(x) && isscalar(x);

defMedFilt = false;
checkMF = @(x) islogical(x) && isscalar(x);

defAllBinFiles = true;
checkABF = checkMF;

defRemArt = false;
checkRA = checkMF;

defTrigWin = 3; % Equivalent to 1.6 ms
checkTW = checkFs;

p.addRequired('dataDir', checkDataDir)
p.addParameter('BinFileName', defOutBinName, checkOutBinName);
p.addParameter('fs', defFs, checkFs);
p.addParameter('ChanNumber', defChan, checkChan);
p.addParameter('MedianFilt', defMedFilt, checkMF);
p.addParameter('AllBinFiles', defAllBinFiles, checkABF);
p.addParameter('RemoveArtifacts', defRemArt, checkRA);
p.addParameter('TriggerWindow', defTrigWin, checkTW);
p.addParameter('verbose', true, checkRA);

parse(p, dataDir, varargin{:});

dataDir = p.Results.dataDir;
outBinName = p.Results.BinFileName;
fs = p.Results.fs;
Nch = p.Results.ChanNumber;
medianFilterFlag = p.Results.MedianFilt;
abfFlag = p.Results.AllBinFiles;
raFlag = p.Results.RemoveArtifacts;
triggSpread = p.Results.TriggerWindow;
verbose = p.Results.verbose;

%% Validation of inputs
iOk = false;
getName = @(x) string(x.name);
getFullName = @(x) string( fullfile(x.folder, x.name) );
fnOpts = {'UniformOutput', false};
searchFileIn = @(inDir, pttrn) dir( fullfile( inDir, pttrn ) );
searchFileHere = @(pttrn) searchFileIn( dataDir, pttrn );
% Search for Recording*.bin files in the given folder
recFiles = searchFileHere('Recording*.bin');
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
if Nrf > 1 && ~abfFlag
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
fNames = arrayfun(getName, recFiles);
rDates = getDates(fNames, 'Recording', '.bin');

% Remove artifacts from triggers - first part: reading the triggers
if raFlag
    trgPttrn = "TriggerSignals*.bin";
    % Search for Trigger*.bin files
    trgFiles = searchFileHere( trgPttrn );
    trgffNames = arrayfun(@(x) getFullName( x ), trgFiles);
    tDates = getDates( trgffNames, 'TriggerSignals', '.bin');
    [tordSubs, ~] = find( rDates == tDates);
    trgFiles = trgFiles(tordSubs);
    if isempty(trgFiles) && verbose
        fprintf(1, 'No trigger files!\n')
        fprintf(1, 'Cannot remove artifacts without the trigger times!\n')
    end
    trigCells = getTriggersFromFiles();
else
    % I feel like this else is going to be needed

end

try
    mem = memory;
    mFlag = true;
catch
    fprintf(1, "Memory pack unavailable!\n")
    mFlag = false;
end
% Overwrite the output file.
fPerm = 'w';
rdf = {@readDataFile, @readDataFileAndMedianFilter};

for cf = 1:Nrf
    fprintf(1, "File: %s\n", fNames(cf));
    % How many bytes is this file?
    fWeight = recFiles(cf).bytes;
    % How many bytes are available in this computer. Occupying 85%
    if mFlag
        mxBytes = 0.8 * mem.MaxPossibleArrayBytes;
    else
        mxBytes = fWeight;
    end
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
        dataBuff = rdf{medianFilterFlag+1}(dfID, Ns, Nch);
        if raFlag
            dataBuff = removeArtifacts(dataBuff, trigCells{cf}, triggSpread);
        end
        ofID = createOrOpenOutputFile(outFullName, fPerm);
        fprintf(1, "Writing...\n")
        if ofID > 2
            fwrite(ofID, dataBuff, "int16")
            fPerm = 'a'; fclose(ofID);
        else
            fprintf(1, 'Error: %d\nCould not write output file!\n', ofID)
            return
        end
    end
    fclose(dfID);
end
iOk = true;
%% Writing useful files
ffoID = fopen(fullfile(dataDir, outBinName + "_fileOrder.txt"), 'w');
fprintf(ffoID, "%s\n", arrayfun(getName, recFiles));
fprintf(ffoID, "%s.bin", outBinName); fclose(ffoID);
save(fullfile(dataDir, outBinName + "_sampling_frequency.mat"), 'fs')

%% Insider functions
    function trigCell = getTriggersFromFiles()
        fIDs = arrayfun(@(x) openDataFile( getFullName(x) ), trgFiles );
        trigs = arrayfun(@(x) readDataFile( x, Inf, 2), fIDs, fnOpts{:} );
        [~] = arrayfun(@(x) fclose(x), fIDs);
        dwObj = cellfun(@(c) arrayfun(@(x) StepWaveform(c(x,:), fs, ...
            'verbose', verbose), 1:size(c, 1) ), trigs, fnOpts{:} );
        trigCell = cellfun(@(x) arrayfun(@(y) y.subTriggers, x, fnOpts{:} ), ...
            dwObj, fnOpts{:} );
    end
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

function [dataBuff] = readDataFile(fID, Ne, Nch)
% As a result to removing the median of the channels, the data will take
% negative values. This is why I need to convert the data to signed 32-bit
% integers before substracting the median. Once the data is median  free,
% then it can be converted back to signed 16-bit integer.
dataBuff = fread(fID, [Nch, Ne], 'uint16=>int32');
dataBuff = int16(dataBuff - median(dataBuff, 2));
end

function fID = createOrOpenOutputFile(fPath,fPerm)
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

function dataBuff = removeArtifacts(dataBuff, triggSubs, triggSpred)
triggSubs = cat(1, triggSubs{:}); triggSubs = sort(triggSubs(:), 'ascend');
triggWin = -triggSpred:triggSpred;
Sw = blackmanharris(numel(triggWin));
for ct = triggSubs'
    for cch = 1:64 % Careful with files with other channel number!
        segm = double(dataBuff(cch, ct + triggWin));
        m = (segm(end) - segm(1))./(triggSpred*2);
        b = segm(1) + m*triggSpred; L = triggWin*m + b;
        segm = int16( round( segm(:) .* ...
            ( ( 1 - Sw(:) ) + Sw(:).*exp(-zscore(double(segm(:))).^2) ) + ...
            Sw(:).*L(:)/2 ) );
        dataBuff(cch, ct + triggWin) = segm;
    end
end
end

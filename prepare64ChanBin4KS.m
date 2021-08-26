function iOk = prepare64ChanBin4KS(dataDir, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Input parsing
p = inputParser;

checkDataDir = @(x) exist(x,'dir');

defOutBinName = 'MedianFiltered_int16';
checkOutBinName = @(x) all([~isempty(x), ischar(x) | isstring(x)]);

p.addRequired('dataDir', checkDataDir)
p.addOptional('BinFileName', defOutBinName, checkOutBinName);

p.parse(dataDir, varargin{:});

dataDir = p.Results.dataDir;
outBinName = p.Results.BinFileName;

%% Validation of inputs
iOk = false;

recFile = dir(fullfile(dataDir, 'Recording*.bin'));
if isempty(recFile)
    fprintf(2, 'No recording file found!\n')
    fprintf(2, 'Make sure you provide a directory with a recording file');
    fprintf(2, ' in it!\n');
    return
end
if numel(recFile) < 2
    % One recording file in directory
    fID = fopen(fullfile(dataDir, recFile.name), 'r');
    if fID < 3
        fprintf(2, 'Unable to open %s!\n', recFile.name)
        fprintf(2, 'Verify that no other process is running with it and')
        fprintf(2, ' try again.\n');
        return
    end
    memObj = memory;
    if 0.85*memObj.MaxPossibleArrayBytes > recFile.bytes
        % It is possible to allocate the whole recording.
        data = fread(fID, Inf, 'uint16'); fclose(fID);
        data = data - 2^15; data = reshape(data, 64, []);
        data = int16(data); expMedian = median(data);
        data = data - expMedian;
        
        fID = fopen(fullfile(dataDir, [outBinName, '.bin']), 'w');
        if fID < 3
            fprintf(2, 'Unable to create output file!\n');
            fprintf(2, 'Please verify the file and directory!\n');
            return
        end % File identifier validation
        fwrite(fID, data, 'int16'); fclose(fID);
        
        fID = fopen(fullfile(dataDir, ['ExperimentMedian.bin']), 'w');
        if fID < 3
            fprintf(2, 'Unable to create median file!\n');
            fprintf(2, 'Please verify the file and directory!\n');
            return
        end % File identifier validation
        fwrite(fID, expMedian, 'int16'); fclose(fID);
        iOk = true;
    else
        % Read a big chunk of data, and write in disk. Use while loop.
        % TODO: Implement this functionality
        
    end % Memory validation
else
    % Deal with more than one recording file
end % More than one Recording* file in the directory

end
function file2bin(fileInput, fileOutput)
% FILEBIN converts a mat, smr, or smrx  file into a binary file to be read
% by Kilosort. It is unclear to me if this step is absolutely necessary,
% but for now it seems helpful. This function takes the path for a data
% file from the lab and saves it in a binary format for Kilosort in the
% (optional) fileOutput path. One limitation of this function is that it
% will not be able to import a analysis.mat file to a binary format. This
% is due a concept incompatibility. The analysis file contains structures
% with different data types and not only signals. This issue includes the
% unability to save the sampling frequency. One solution to this problem
% might be that the sampling frequency is at the very beguinning of the
% file, in the first channel, in the first position.
    
if isstring(fileInput)
    fileInput = char(fileInput);
end
fprintf('Reading %s file...\n',fileInput)
if exist(fileInput,'file')
    [flDir, baseName, flEx] = fileparts(fileInput);
    if nargin < 2
        fileOutput = fullfile(flDir,[baseName,'.bin']);
    end
    switch flEx
        % Depending on the file extension, the function will import and
        % save the binary file in different manners
        case {'.smr','.smrx'}
            fprintf('spike2 file format recognized\n')
            fID = fopen(fullfile(flDir,[baseName,flEx]),'r');
            SONXimport(fID,'bin');
        case '.mat'
            fprintf('.mat file recognized\n')
            if contains(baseName,'analysis')
                fprintf('The analysis file is not meant to be imported\n')
                baseName = strsplit(baseName,'analysis');
                baseName = baseName{1};
                mat2bin(fullfile(flDir,[baseName,'.mat']),fileOutput)
            end
            fprintf('Importing %s...\n',[baseName,flEx])
        otherwise
            fprintf('This input file is not recognized.')
            return
    end
else
    fprintf('The input file was not found')
end
end

function mat2bin(fileFullPath, fileOutput)
if ~exist(fileOutput,'file')
    [foDir,foBaseName,foExt] = fileparts(fileOutput);
    dataLoader = UMSDataLoader(fileFullPath);
    data = dataLoader.getDataMatrix;
    fs = dataLoader.SamplingFrequency;
    fsRow = zeros(1,dataLoader.Nch);
    fsRow(1) = fs;
    data = [fsRow;data]';
    try
        if ~strcmp(foExt,'.bin')
            fileOutput = fullfile(foDir,[foBaseName,'.bin']);
        end
        fid_out = fopen(fileOutput,'w');
        fwrite(fid_out,data,'*int16');
        fclose(fid_out);
    catch
        disp('There was an error creating/writing the file')
    end
else
    disp('The file exists already.') 
end
end
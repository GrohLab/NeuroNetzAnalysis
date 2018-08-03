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
if nargin == 2
    % The file output was given.
    fprintf('Reading %s file...\n',fileInput)
    if exist(fileInput,'file')
        [flDir, baseName, flEx] = fileparts(fileInput);
        switch flEx
            % Depending on the file extension, the function will import and
            % save the binary file in different manners
            case {'.smr','.smrx'}
                fprintf('spike2 file format recognized')
                crrDir = pwd;
                importSMR([baseName,flEx],flDir,1);
                cd(crrDir)
                mat2bin(fullfile(flDir,[baseName,'.mat']),fileOutput)
            case '.mat'
                fprintf('.mat file recognized')
                if contains(baseName,'analysis')
                    fprintf('The analysis file is not meant to be imported\n')
                    baseName = strsplit(baseName,'analysis');
                    baseName = baseName{1};
                end
                fprintf('Importing %s...\n',[baseName,flEx])
                
            otherwise
                fprintf('This input file is not recognized.')
        end
    else
        fprintf('The input file was not found')
    end
elseif nargin < 2
    
end

end

function mat2bin(fileFullPath, fileOutput)
if ~exist(fileOutput,'file')
    [foDir,foBaseName,foExt] = fileparts(fileOutput);
    [fiDir,fiBaseName,~] = fileparts(fileFullPath);
    dataLoader = UMSDataLoader(fileFullPath);
    data = dataLoader.getDataMatrix;
    fs = dataLoader.SamplingFrequency;
    fsRow = zeros(1,dataLoader.Nch);
    fsRow(1) = fs;
    data = [fsRow;data]';
    try
    fid_out = fopen(fileOutput,'w');
    fwrite(fid_out,data,'double')
    fclose(fid_out);
    catch
        disp('There was an error creating/writing the file')
    end
else
    ovwrt = input('The file exists already. Do you wish to overwrite it?','s');
    if 
end
end
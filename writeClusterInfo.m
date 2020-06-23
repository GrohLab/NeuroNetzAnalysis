function iOk = writeClusterInfo(clInfo, fileName, owFlag)
iOk = false;
askOw = true;
if exist('owFlag','var') && logical(owFlag)
    askOw = false;
end

% If the file exists, ask if it should be overwritten.
if exist(fileName,'file') && askOw
    owAns = questdlg(sprintf('%s exists. Overwirte?',fileName),...
        'Overwrite','Yes','No','No');
    if strcmpi(owAns,'no')
        fprintf(1, 'No file writen.\n');
        return
    end
end
% Table given
if ~istable(clInfo)
    fprintf(1,'The given cluster table is not a MATLAB table. ')
    fprintf(1,'Cannot continue.\n')
    return
end

% Validation of the input file
[dataDir, ~, ~] = fileparts(fileName);
if isempty(dataDir)
    dataDir = pwd;
    fprintf(1, 'Warning: writing file in %s', dataDir);
end

fID = fopen(fileName, 'w');
% Write header
varNames = string(clInfo.Properties.VariableNames);
for cvn = 1:length(varNames)-1
    fprintf(fID, '%s\t', varNames(cvn));
end
fprintf(fID, '%s\n', varNames(cvn+1));
% Write data
[Ncl, Nv] = size(clInfo);
for ccl = 1:Ncl
    sep = '\t';
    for cv = 1:Nv
        if cv == Nv
            sep = '\n';
        end
        switch varNames(cv)
            case {"id","KSLabel","group","NeuronType"}
                fprintf(fID,['%s',sep],clInfo{ccl,cv}{1});
            case "firing_rate"
                fprintf(fID,['%f spk/s',sep],clInfo{ccl,cv});
            otherwise
                fprintf(fID,['%f',sep],clInfo{ccl,cv});
        end
    end
end
fclose(fID);
iOk = true;
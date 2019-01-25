function [iok] =...
    savePopConfigFile(conFiNa, configStruct)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if exist(conFiNa,'file')
    warning('The file exists!')
    rwt = questdlg(sprintf('Rewrite? (%s)',conFiNa),'File conflict','Yes',...
        'No','No');
    if strcmp(rwt,'No')
        iok = 1;
        fprintf('File not written.\n')
        return;
    end
end
% Write the POPULATION configuration file.
fID = fopen(conFiNa,'w'); % Overwrite any existing file
% Create a "header" (only date and time)
fprintf(fID,'date:\t%s',datetime('now'));
writeNSignalNames('cellType(s):',fID,configStruct,'CellType')
% Viewing windows
fprintf(fID,'viewWind:\t%f\n',configStruct.ViewWindow);
% Binning time
fprintf(fID,'binSz:\t%f\n',configStruct.BinSize);
% Trigger name and considered edge
fprintf(fID,'t:\t%s\t%d',configStruct.Trigger.Name,...
    configStruct.Trigger.Edge);
% Exclude signals
writeNSignalNames('e:',fID,configStruct,'Exclude')
% Ignore signals
writeNSignalNames('i:',fID,configStruct,'Ignore')
% Conditioning windows
fprintf(fID,'w:');
if isempty(configStruct.ConditionWindow)
    fprintf(fID,'\tnone\t%d %d\n',0,0);
else
    for cws = 1:length(configStruct.ConditionWindow.Names)
        fprintf(fID,'\t%s\t%f',configStruct.ConditionWindow.Names{cws},...
            configStruct.ConditionWindow.ConditionWindow);
    end
    fprintf(fID,'\n');
end
iok = fclose(fID);
if ~iok
    fprintf('File successfully written!\n')
    fprintf('%s',conFiNa)
else
    fprintf('Something went wrong writing the file...\n')
    fprintf('It is possible that it contains no data...\n')
end
end

function writeNSignalNames(preamb, fID, cStr, fiNa)
fprintf(fID,preamb);
for cs = 1:length(cStr.(fiNa))
    fprintf(fID,'\t%s',cStr.(fiNa){cs});
end
fprintf(fID,'\n');
end


% fprintf(fID,'cellType(s):');
% for cct = 1:length(configStruct.CellType)
%     fprintf(fID,'\t%s',configStruct.CellType{cct});
% end
% fprintf(fID,'\n');
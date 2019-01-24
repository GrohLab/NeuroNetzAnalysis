function [iok] =...
    savePopConfigFile(conFiNa, configStruct)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if exist(conFiNa,'file')
    warning('The file exists!')
    rwt = inputdlg(sprintf('Rewrite? (%s)',conFiNa),'File conflict');
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
fprintf(fID,'cellType(s):');
for cct = 1:length(configStruct.CellType)
    fprintf(fID,'\t%s',configStruct.CellType{cct});
end
fprintf(fID,'\n');
% Viewing windows
fprintf(fID,'viewWind:\t%f\n',configStruct.ViewWindow);
% Binning time
fprintf(fID,'binSz:\t%f\n',configStruct.BinSize);
% Trigger name and considered edge
fprintf(fID,'t:\t%s\t%d',configStruct.Trigger.Name,...
    configStruct.Trigger.Edge);
% Exclude signals
fprintf(fID,'e:');
for ces = 1:length(configStruct.Exclude)
    fprintf(fID,'\t%s',configStruct.Exclude{ces});
end
fprintf(fID,'\n');
% Ignore signals
fprintf(fID,'e:');
for cis = 1:length(configStruct.Ignore)
    fprintf(fID,'\t%s',configStruct.Ignore{cis});
end
fprintf(fID,'\n');
end
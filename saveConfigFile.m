function [iok] =...
    saveConfigFile(conFiNa, configStruct)
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
% Create a header
fprintf(fID,'date:\t%s',datetime('now'));
fprintf(fID,'cellType(s):');
for cct = 1:length(configStruct.CellType)
    fprintf(fID,'\t%s',configStruct.CellType{cct});
end
fprintf(fID,'\n');
fprintf(fID,'t:\t%s\t%d','yes',1);
end


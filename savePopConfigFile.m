function [iok] =...
    savePopConfigFile(conFiNa, configStruct)
%SAVEPOPCONFIGFILE saves the configuration file for a population analysis.
%   The configuration structure in _configStruct_ is created at the
%   beguinning of the analysis if there's a new configuration. There should
%   be a configuration history to easily select a previous configuration.
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
fprintf(fID,'date:\t%s\n',datetime('now'));
writeNSignalNames('cellType(s):',fID,configStruct,'CellType')
% Viewing windows
fprintf(fID,'viewWind:\t%f %f\n',configStruct.ViewWindow(1),...
    configStruct.ViewWindow(2));
% Binning time
fprintf(fID,'binSz:\t%f\n',configStruct.BinSize);
% Trigger name and considered edge
fprintf(fID,'t:\t%s %d\n',configStruct.Trigger.Name,...
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
    for cws = 1:length(configStruct.ConditionWindow)
        fprintf(fID,'\t%s\t%f %f',configStruct.ConditionWindow(cws).Name,...
            configStruct.ConditionWindow(cws).Window(cws,1),...
            configStruct.ConditionWindow(cws).Window(cws,2));
    end
    fprintf(fID,'\n');
end
iok = fclose(fID);
if ~iok
    fprintf('File successfully written!\n')
    fprintf('%s\n',conFiNa)
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
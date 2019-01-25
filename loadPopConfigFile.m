function [configStruct] = loadPopConfigFile(conFiNa)
%LOADPOPCONFIGFILE gets the information in the text file and returns a
%configuration structure in the same manner as the savePopConfigFile
%accepts it.
%   The output structure contains all the information necessary to perform
%   a population analysis. The validation of the file ending should be
%   performed before entering the 
if exist(conFiNa,'file')
    fID = fopen(conFiNa,'r');
    while ~feof(fID)
        ln = fgetl(fID);
        tsl = strsplit(ln,'\t');
        switch tsl{1}
            case 'date:'
                dt = strsplit(tsl{2},' ');
                fprintf('Analysis performed the %s at %s\n',dt{1},dt{2})
            case 'cellType(s):'
                Nct = length(tsl);
                configStruct.CellType = tsl(2:Nct)';
            case 'viewWind:'
                timeLapse = str2double(strsplit(tsl{2},' '));
                configStruct.ViewWindow = timeLapse;
            case 'binSz:'
                configStruct.BinSize = str2double(tsl{2});
            case 't:'
                trigInfos = strsplit(tsl{2},' ');
                triggerStruct.Name = trigInfos{1};
                triggerStruct.Edge = str2double(trigInfos{2}) == 1;
                configStruct.Trigger = triggerStruct;
            case 'e:'
                Ne = length(tsl);
                configStruct.Exclude = tsl(2:Ne)';
            case 'i:'
                Ni = length(tsl);
                configStruct.Ignore = tsl(2:Ni)';
            case 'w:'
                Nw = (length(tsl) - 1)/2;
                for cw = 1:Nw
                    auxStruct.Name = tsl{cw*2};
                    cwd = str2double(strsplit(tsl{cw*2 + 1}));
                    auxStruct.Window = cwd;
                    configStruct.ConditionWindow(cw) = auxStruct;
                end
            otherwise
                fprintf('New configuration parameter found!\n')
                fprintf('Please program this function.\n')
        end
    end
else
    fprintf('File not found.\n')
    configStruct = [];
end
end
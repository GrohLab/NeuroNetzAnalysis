function [rollerposition, tTimes] = readRollerPositionsFile(filepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tTimes = [];
fID = fopen(filepath, 'r');
rollerline = textscan(fID, '%d,2021-%*d-%*dT%d:%d:%f+02:00');
fclose(fID);
if all(cellfun(@(x) ~isempty(x), rollerline))
    % Roller file corresponding to the Bonsai timestamps
    rollerx = int16(rollerline{1}); rollert = string(rollerline{2})+":"+...
        rollerline{3}+":"+rollerline{4}; rollert = duration(rollert);
    rollert.Format = 'hh:mm:ss.SSSSSSS';
    rollerposition = table(rollerx, rollert, ...
        'VariableNames', {'RollerX','RollerT'});
else
    fnOpts = {"UniformOutput", false};
    ops = delimitedTextImportOptions("NumVariables", 2);
    ops.DataLines = [1, Inf];
    ops.Delimiter = ",";
    ops.VariableNames = ["RoT", "Time_mus"];
    ops.VariableTypes = ["string", "double"];
    ops.ExtraColumnsRule = "ignore";
    ops.EmptyLineRule = "read";
    % Roller file corresponding to Arduino's micros() function.
    rollTrigTimes = readtable(filepath, ops);
    rollTrigTimes.Time_mus = unwrap(rollTrigTimes.Time_mus, 2^32);
    % Is there any 
    trigFlag = isnan(str2double(rollTrigTimes.RoT));
    trigLetter = unique(rollTrigTimes.RoT(trigFlag)); 
    if ~isempty(trigLetter)
        trigNum = size(trigLetter, 1);
        fprintf(1, "Found %d trigger(s): %s\n", trigNum,...
            sprintf("%s ",trigLetter))
        % Comparing the content of the table with all found letters
        % corresponding to specific triggers
        tFlags = arrayfun(@(x) ismember(rollTrigTimes.RoT(trigFlag), x),...
            trigLetter', fnOpts{:}); 
        % Reading those times which correspond only to the specific trigger
        tTimes = cellfun(@(x) rollTrigTimes.Time_mus(x), tFlags, fnOpts{:});
    end
end
end
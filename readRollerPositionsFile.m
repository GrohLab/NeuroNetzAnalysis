function [rollerposition, tTimes, rollTrigTimes] = readRollerPositionsFile(filepath)
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
    rollTrigTimes.Time_mus = unwrap(rollTrigTimes.Time_mus, 2^31);
    % Searching for well-written trigger interruptions
    trigFlag = isnan(str2double(rollTrigTimes.RoT));
    trigLetter = unique(rollTrigTimes.RoT(trigFlag));
    % Rescueing the lost trigger interruptions
    nanSubs = find(isnan(rollTrigTimes.Time_mus));
    fID = fopen(filepath, "r"); ln = 1;
    for cns = nanSubs'
        while ln < cns
            fgetl(fID); ln = ln + 1;
        end
        strLn = fgetl(fID); ln = ln + 1;
        strCell = strsplit(strLn, ","); chrepFlag = strCell{2} > 57;
        refTm = num2str(mean(rollTrigTimes.Time_mus(cns+[1;-1])));
        rollTrigTimes.RoT(cns) = string(strCell{2}(chrepFlag));
        trigFlag(cns) = true;
        strCell{2}(chrepFlag) = refTm(chrepFlag);
        rollTrigTimes.Time_mus(cns) = str2double(strCell{2});
    end
    [~] = fclose(fID);
    % Output for the roller positions
    rollerposition = zeros(sum(~trigFlag), 2);
    rollerposition(:,2) = rollTrigTimes.Time_mus(~trigFlag);
    rollerposition(:,1) = unwrap(str2double(rollTrigTimes.RoT(~trigFlag)),...
        2^15);
    % Output for the trigger times as measured by the rotary decoder
    if ~isempty(trigLetter)
        trigNum = size(trigLetter, 1);
        fprintf(1, "Found %d trigger(s): %s\n", trigNum,...
            sprintf("%s ",trigLetter))
        % Comparing the content of the table with all found letters
        % corresponding to specific triggers
        tFlags = arrayfun(@(x) ismember(rollTrigTimes.RoT(trigFlag), x),...
            trigLetter', fnOpts{:}); trigSubs = find(trigFlag);
        % Reading those times which correspond only to the specific trigger
        tTimes = cellfun(@(x) rollTrigTimes.Time_mus(trigSubs(x)), tFlags,...
            fnOpts{:});
    end
end
end
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
    ops.VariableNames = ["RoT", "Time_us"];
    ops.VariableTypes = ["string", "double"];
    ops.ExtraColumnsRule = "ignore";
    ops.EmptyLineRule = "read";
    % Roller file corresponding to Arduino's micros() function.
    rollTrigTimes = readtable(filepath, ops);
    % Accounting for long int format in the Arduino clock.
    rollTrigTimes.Time_us = unwrap(rollTrigTimes.Time_us, 2^31);
    % Searching for well-written trigger interruptions
    trigFlag = any([isnan(str2double(rollTrigTimes.RoT)),...
        isnan(rollTrigTimes.Time_us)], 2);
    trigLetter = unique(rollTrigTimes.RoT(trigFlag));
    % Rescuing ill-written trigger time interruptions
    nanSubs = find(isnan(rollTrigTimes.Time_us));
    fID = fopen(filepath, "r"); ln = 1;
    for cns = nanSubs'
        while ln < cns
            prevLn = fgetl(fID); ln = ln + 1; %#ok<NASGU>
        end
        strLn = fgetl(fID); ln = ln + 1;
        if ~isempty(strLn)
            strCell = strsplit(strLn, ",");
            if size(strCell,2) == 2
                % Character outside numeric ASCII range.
                chrepFlag = strCell{2} > 57;
                % Taking the mean from 6 cells around the found trigger
                % time.
                refTm = num2str(round(mean(rollTrigTimes.Time_us(cns+(-3:3)'),...
                    'omitnan'))); 
                % Writting the found trigger character and logic flag.
                rollTrigTimes.RoT(cns) =...
                    string(strCell{2}(chrepFlag)); trigFlag(cns) = true;
                % "Mistake" correction cases:
                if length(refTm) == length(strCell{2})
                    % If the letter is in the middle of the microseconds,
                    % replace the character with the mean value from all 6
                    % previous cells
                    strCell{2}(chrepFlag) = refTm(chrepFlag);
                elseif length(strCell{2}) > length(refTm)
                    % If the 
                    strCell{2}(chrepFlag) = [];
                else
                    strCell{2} = refTm;
                end
                rollTrigTimes.Time_us(cns) = str2double(strCell{2});
            else
                % There are more commas in the line. Something like an
                % interrupted interruption, new line in the time string,
                % and continued writing the time string.                
                fprintf(1, 'Need user interaction!\n');
            end
        else
            % Empty line--maybe double enter from interrupted interruption
            rollTrigTimes(cns,:) = [];
            trigFlag(cns) = [];
        end
    end
    [~] = fclose(fID);
    % Output for the roller positions
    rollerposition = zeros(sum(~trigFlag), 2);
    rollerposition(:,2) = rollTrigTimes.Time_us(~trigFlag);
    rollerposition(:,1) = unwrap(str2double(rollTrigTimes.RoT(~trigFlag)),...
        2^15);
    % Output for the trigger times as measured by the rotary decoder
    if ~isempty(trigLetter)
        % Searching for ill-written trigger ID interruptions
        lettFlag = strlength(trigLetter) == 1;
        trigLetter(~lettFlag) = []; trigNum = size(trigLetter, 1);
        fprintf(1, "Found %d trigger(s): %s\n", trigNum,...
            sprintf("%s ",trigLetter))
        % Comparing the content of the table with all found letters
        % corresponding to specific triggers
        tFlags = arrayfun(@(x) ismember(rollTrigTimes.RoT(trigFlag), x),...
            trigLetter', fnOpts{:}); trigSubs = find(trigFlag);
        % Reading those times which correspond only to the specific trigger
        tTimes = cellfun(@(x) rollTrigTimes.Time_us(trigSubs(x)), tFlags,...
            fnOpts{:});
    end
end
end
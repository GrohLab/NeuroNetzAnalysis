function [rollerposition, tTimes, rollTrigTimes] = readRollerPositionsFile(filepath)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tTimes = [];
fID = fopen(filepath, 'r');
rollerline = textscan(fID, '%d,2021-%*d-%*dT%d:%d:%f+02:00');
fclose(fID);

%% Auxiliary functions
refTime = 0; refRoT = 'C'; nxtRoT = '';nxtTime = 0;

    function placeTriggerInNextRow()
        % Save the time fragment from the position, if any.
        fetchNextRoT()
        if ismissing(nxtRoT)
            nxtRoT = [];
        end
        % Setting the trigger flag to false, as the position is not lost,
        % and the next line true.
        trigFlag(cns) = false; trigFlag(cns+1) = true;
        rollTrigTimes.RoT(cns+1) = strCell{clFlag}(chRepFlag);
        rollTrigTimes.Time_us(cns+1) = str2double(strCell{end});
        % Fixing the current line by concatenating the time stamp of the
        % roller position
        strCell{clFlag}(chRepFlag) = [];
        strCell{clFlag} = cat(2, strCell{clFlag}, char(nxtRoT));
        % rollTrigTimes.RoT(cns) = string(strCell{1});
        rollTrigTimes.Time_us(cns) = str2double(strCell{clFlag});
    end

    function placeTriggerInCurrentRow()
        fetchNextRoT()
        if ismissing(nxtRoT)
            nxtRoT = [];
        end
        % Setting the trigger flag to false, as the position is not lost,
        % and the next line true.
        trigFlag(cns) = true; trigFlag(cns+1) = false;
        % Write the trigger letter and time in the following line.
        % rollTrigTimes(cns+1,:)
        if any(strCell{2} == 45)
            % The position is negative. Easier to split.
            tmPos = strsplit(strCell{2}, '-');
            tmPos{1} = cat(2, tmPos{1}, char(nxtRoT));
            rollTrigTimes.RoT(cns+1) = cat(2, '-', tmPos{2});
            rollTrigTimes.Time_us(cns) = str2double(tmPos{1});
        else
            % Difficult to separate position and time. In development
            refRoT = createRefRoT();
            rollTrigTimes.Time_us(cns) = strCell{end};
            rollTrigTimes.RoT(cns+1) = num2str(refRoT);
            % TODO: Implement a way to recognise the time from the position
            % numPatDm = (strCell{2} - num2str(refRoT)') == 0;
        end
        rollTrigTimes.Time_us(cns+1) = str2double(strCell{end});
    end

    function correctTimeAndReplacePosWithTrigger()
        createRefTime(); refTime = char(string(round(refTime)));
        rollTrigTimes.RoT(cns) = strCell{clFlag}(chRepFlag);
        if length(refTime) == length(strCell{clFlag})
            % If the letter is in the middle of the microseconds, replace
            % the character with the mean value from all 20 cells
            strCell{clFlag}(chRepFlag) = refTime(chRepFlag);
        elseif length(strCell{clFlag}) > length(refTime)
            % This might mean that the trigger letter is at the end of the
            % line. The solution is just to delete the trigger letter.
            strCell{clFlag}(chRepFlag) = [];
        else
            % Strange case. I cannot imagine what could have happened here.
            strCell{clFlag} = refTime;
        end
        rollTrigTimes.Time_us(cns) = str2double(strCell{clFlag});
        % rollTrigTimes(cns,:)
        % trigFlag(cns)
    end

    function correctPosAndPlaceTrigNextRow()
        % Trigger is in the position place and need to validate more.
        trigStrX = find(chRepFlag);
        if trigStrX == size(chRepFlag,2)
            % Trigger at the end of the string. Need to search in the next
            % line for correction clues.
            fetchNextRoT()
            createRefRoT()
            if ismissing(nxtRoT) ||...
                    strlength(nxtRoT) < strlength(string(refRoT))
                trigFlag(cns) = false; trigFlag(cns+1) = true;
                % Bring up the original position time
                rollTrigTimes.Time_us(cns) =...
                    rollTrigTimes.Time_us(cns+1);
                % The position writing was interrupted by the trigger
                rollTrigTimes.RoT(cns+1) = strCell{clFlag}(chRepFlag);
                rollTrigTimes.Time_us(cns+1) = str2double(strCell{end});
                % Fixing the current line by concatenating the time stamp
                % of the roller position
                strCell{clFlag}(chRepFlag) = [];
                strCell{clFlag} = cat(2, strCell{clFlag}, char(nxtRoT));
                rollTrigTimes.RoT(cns) = strCell{clFlag};
            else
                % Trigger char was written at the end of the position, but
                % without it's own time.
                rollTrigTimes.RoT(cns) = strCell{clFlag}(chRepFlag);
            end
        else
            % Trigger character in the middle of the position but without
            % it's own time.
            rollTrigTimes.RoT(cns) = strCell{clFlag}(chRepFlag);
        end
    end

    function fetchNextRoT()
        nxtRoT = rollTrigTimes.RoT(cns+1);
    end

    function fetchNextTime()
        nxtTime = rollTrigTimes.Time_us(cns+1);
    end

    function refRoT = createRefRoT()
        % Create reference roller position number.
        refRoT = mean(str2double(rollTrigTimes.RoT(cns+([(-10:-1)';...
            (1:10)']))),'omitnan');
    end

    function createRefTime()
        % Create reference time.
        refTime = mean(rollTrigTimes.Time_us(cns+([(-10:-1)';(1:10)'])),...
            'omitnan');
    end

%% Reading file
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
    ops.MissingRule = "fill";
    ops.VariableNames = ["RoT", "Time_us"];
    ops.VariableTypes = ["string", "double"];
    ops.ExtraColumnsRule = "ignore";
    ops.EmptyLineRule = "read";
    ops.ConsecutiveDelimitersRule = "join";
    % Roller file corresponding to Arduino's micros() function.
    rollTrigTimes = readtable(filepath, ops);
    % Searching for well-written trigger interruptions
    trigFlag = any([isnan(str2double(rollTrigTimes.RoT)),...
        isnan(rollTrigTimes.Time_us)], 2);
    trigLetter = unique(rollTrigTimes.RoT(trigFlag));
    % Rescuing ill-written rows
    nanSubs = find(trigFlag & strlength(rollTrigTimes.RoT) > 1);
    fID = fopen(filepath, "r"); ln = 1;
    for cns = nanSubs'
        while ln < cns
            prevLn = fgetl(fID); ln = ln + 1; %#ok<NASGU>
        end
        strLn = fgetl(fID); ln = ln + 1;
        % Split the line by commas.
        strCell = strsplit(strLn, ",");
        Ncell = size(strCell,2);
        % Searching for the cell which has a letter.
        clFlag = cellfun(@(x) isnan(str2double(x)), strCell);
        % Character outside numeric ASCII range on the second cell.
        if any(clFlag)
            chRepFlag = strCell{clFlag} > 57;
        end
        clX = find(clFlag);
        if Ncell > 3
            % Not considered so far. Better prompt the user to check the
            % file!
            fprintf(1, 'Check the considered file!\n')
            fprintf(1, 'Line %d:%s\n', cns, strLn);
            return;
        elseif Ncell > 2
            % Line contains 2 commas separating 3 strings.
            if clX > 1
                % Found trigger char in the second string. The position
                % string is possibly cut.
                placeTriggerInNextRow()
            else
                % Found in the position string.
                placeTriggerInCurrentRow()
            end
        elseif Ncell > 1
            % Line contains 1 comma separating 2 strings.
            if clX > 1
                % The character is in the second cell.
                correctTimeAndReplacePosWithTrigger()
            else
                % The character is in the first cell.
                correctPosAndPlaceTrigNextRow()
            end
        else
            % Either one string or an empty line. Checking position/trigger
            % and time. Very unlikely to enter in any of these ifs
            if ismissing(rollTrigTimes.RoT(cns))
                createRefRoT()
                rollTrigTimes.RoT(cns) = string(refRoT);
            end
            if isnan(rollTrigTimes.Time_us(cns))
                createRefTime()
                rollTrigTimes.Time_us(cns) = refTime;
            end
        end
    end
    [~] = fclose(fID);
    % Output for the roller positions
    rollerposition = zeros(sum(~trigFlag), 2);
    % Accounting for long int format in the Arduino clock.
    rollTrigTimes.Time_us = unwrap(rollTrigTimes.Time_us, 2^31);
    rollerposition(:,1) = unwrap(str2double(rollTrigTimes.RoT(~trigFlag)),...
        2^15);
    rollerposition(:,2) = rollTrigTimes.Time_us(~trigFlag);
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

%{
if size(strCell,2) > 1
                
                if size(strCell,2) < 3
                    % Writting the found trigger character and logic flag.
                    rollTrigTimes.RoT(cns) =...
                        string(strCell{clFlag}(chRepFlag));
                    trigFlag(cns) = true;
                    % If the trigger letter was written on the time string
                    if find(clFlag) == 2
                        % Taking the mean from 6 cells around the found
                        % trigger time. rollTrigTimes(cns+(-3:3)',:)
                        refNum = num2str(round(mean(rollTrigTimes{...
                            cns+(-3:3)',clFlag}, 'omitnan')));
                        % "Mistake" correction cases:
                        if length(refNum) == length(strCell{clFlag})
                            % If the letter is in the middle of the
                            % microseconds, replace the character with the
                            % mean value from all 6 previous cells
                            strCell{clFlag}(chRepFlag) = refNum(chRepFlag);
                        elseif length(strCell{clFlag}) > length(refNum)
                            % This might mean that the trigger letter is at
                            % the end of the line. The solution is just to
                            % delete the trigger letter.
                            strCell{clFlag}(chRepFlag) = [];
                        else
                            % Strange case. I cannot imagine what could
                            % have happened here.
                            strCell{clFlag} = refNum;
                        end
                    else
                        % The trigger letter was written on the position
                        % part of the line. There's nothing more to be done
                        % here. Replaced the position time stamp with the
                        % trigger.
                        continue
                    end
                else
                    % There are more commas in the line. Something like an
                    % interrupted interruption, new line in the time
                    % string, and continued writing the time string.
                    
                    % Save the time fragment from the position, if any.
                    nxtRoT = rollTrigTimes.RoT(cns+1);
                    if ismissing(nxtRoT)
                        nxtRoT = [];
                    end
                    % Setting the trigger flag to false, as the
                    % position is not lost, and the next line true.
                    trigFlag(cns) = false; trigFlag(cns+1) = true;
                    % Write the trigger letter and time in the
                    % following line. rollTrigTimes(cns+1,:)
                    rollTrigTimes.RoT(cns+1) = strCell{clFlag}(chRepFlag);
                    rollTrigTimes.Time_us(cns+1) = str2double(strCell{3});
                    % Fixing the current line by concatenating the time
                    % stamp of the roller position
                    strCell{clFlag}(chRepFlag) = [];
                    strCell{clFlag} = cat(2, strCell{clFlag}, char(nxtRoT));
                end
                switch find(clFlag)
                    case 1
                        rollTrigTimes.RoT(cns) = strCell{clFlag};
                    case 2
                        rollTrigTimes.Time_us(cns) =...
                            str2double(strCell{clFlag});
                    otherwise
                        fprintf(1, 'A case worth seeing!\n')
                end
            else
                fprintf(1, 'When we only have a time fragment\n');
            end
%}
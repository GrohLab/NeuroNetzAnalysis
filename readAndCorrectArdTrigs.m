function [Conditions_arduino] = readAndCorrectArdTrigs(dataDir)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%% Auxiliary variables and functions
fnOpts = {"UniformOutput", false};
dateFormStr = 'yyyy-MM-dd''T''HH_mm_ss';
% Function to extract the date from the file name. 
% Extract machine readable 19-December-2021 15:32:12 from
% Roller_position2021-12-19T15_32_12.csv, for example. "Roller_position"
% would be the baseName, and an array with structures from dir function.
    function fDates = getDates(fileNames, baseName)
        getDate = @(x, y) extractBefore(extractAfter(x, y), ".");
        fDates = arrayfun(@(x) getDate(x.name, baseName), fileNames,...
            fnOpts{:});
        fDates = cellfun(@(x) datetime(x, 'InputFormat', dateFormStr),...
            fDates);
    end
% Function to get the edges from the trigger signals in the trigger file.
    function tSubs = getSubsFromTriggers(trig)
        % Unstable line! Display function not working!
        tObj = arrayfun(@(x) StepWaveform(trig(x,:), 3e4, units='',...
            title='', verbose=false), [1;2]);
        tSubs = arrayfun(@(x) x.subTriggers, tObj, fnOpts{:});
        tObj(2).MinIEI = 1; fot = tObj(2).FirstOfTrain;
        % if the laser signal has no frequency, fot will be populated by
        % zeroes.
        if sum(fot)
            tSubs(2) = {tSubs{2}(fot,:)};
        end
    end
%% Validation section
% Given folder exists?
Conditions_arduino = struct();
if ~exist(dataDir, "dir")
    fprintf(1, "Given directory doesn''t exist!\n")
    return
end
% Search pattern for the rest of the files
sp = fullfile(dataDir, "**");
% Does the folder have files with roller positions?
rpFiles = dir(fullfile(sp, "Roller_position*.csv"));
if isempty(rpFiles)
    fprintf(1, "No roller position files found in:\n%s\n", dataDir)
    return
end
% Does the experiment contain intan trigger signals?
itFiles = dir(fullfile(sp, "TriggerSignals*.bin"));
if isempty(itFiles)
    fprintf(1, "No RHD Intan trigger files found in:\n%s\n", dataDir)
    return
end
% Are the files corresponding to each other?
% Checking for number of files
if numel(rpFiles) ~= numel(itFiles)
    fprintf(1, "Discrepancy between files!\n")
    fprintf(1, "# roller files: %d, # trigger files: %d\n",...
        numel(rpFiles), numel(itFiles));
    return
end
% Extract the date from the file names. Do the names correspond?
rpDates = getDates(rpFiles, "Roller_position");
itDates = getDates(itFiles, "TriggerSignals");
if ~all(rpDates - itDates)
    fprintf(1, "The file names do not correspond!\n")
    fprintf(1, "Please check and maybe even organise the files in this ")
    fprintf(1, "folder\n"); fprintf(1, "Roller files:\n");
    fprintf(1, "%s\n", arrayfun(@(x) string(x.name), rpFiles))
    fprintf(1, "Trigger files:\n");
    fprintf(1, "%s\n", arrayfun(@(x) string(x.name), itFiles))
    return
end

%% Searching for arduino instabilities
for cfp = 1:numel(rpFiles)
    % Read the roller position file
    [rp, tTimes] = readRollerPositionsFile(...
        fullfile(rpFiles(cfp).folder, rpFiles(cfp).name));
    % Converting the trigger times from microseconds to seconds
    tTimes(1,:) = cellfun(@(x) (x - rp(1,2))/1e6, tTimes(1,:), fnOpts{:});
    % Deleting noise entries to the arduino (times with less than 1 second
    % difference)
    delFlags = cellfun(@(x) [false;diff(x) < 1], tTimes(1,:), fnOpts{:});
    tTimes(1,:) = cellfun(@(x,y) x(~y), tTimes(1,:), delFlags, fnOpts{:});
    % Arduino times
    atTimes = tTimes(1,:);
    % Read trigger signals
    tfName = fullfile(itFiles(cfp).folder,itFiles(cfp).name);
    trig = readTriggerFile(tfName);
    % Subs
    tSubs = getSubsFromTriggers(trig);
    itTimes = cellfun(@(x) x./3e4, tSubs, fnOpts{:});
   
end
end

function trig = readTriggerFile(tfName)
fID = fopen(tfName,'r'); trig = fread(fID, [2,Inf], 'uint16=>int32');
[~] = fclose(fID); trig = int16(trig - median(trig, 2));
trig(1,:) = -trig(1,:);
end


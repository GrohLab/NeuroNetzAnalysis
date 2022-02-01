function [atTimes, itTimes, Conditions_arduino] = readAndCorrectArdTrigs(dataDir)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%% Auxiliary variables and functions
fnOpts = {"UniformOutput", false};
dateFormStr = 'yyyy-MM-dd''T''HH_mm_ss';
derv = @(x) diff(x(:,1));
dsmt = @(x,y) distmatrix(x, y, 2);
erOpts = {"ErrorHandler", @falseLaserDetection};
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
% Function to get arduino trigger times
    function [atTimes, atNames] = getArduinoTriggers(rpFilePath)
        % Read the roller position file
        [rp, tTimes] = readRollerPositionsFile(rpFilePath);
        % Converting the trigger times from microseconds to seconds
        tTimes(1,:) = cellfun(@(x) (x - rp(1,2))/1e6, tTimes(1,:), fnOpts{:});
        % Deleting noise entries to the arduino (times with less than 1 second
        % difference)
        delFlags = cellfun(@(x) [false;diff(x) < 1], tTimes(1,:), fnOpts{:});
        tTimes(1,:) = cellfun(@(x,y) x(~y), tTimes(1,:), delFlags, fnOpts{:});
        % Arduino times
        atTimes(cfp,:) = tTimes(1,:); atNames = string(tTimes(2,:));
    end
% Function to compute 2 matrices: similarity between trigger intervals
% (ddm) and similarity between trigger times (dm)
    function [ddm, dm] = computeSimilarityMatrices()
        atDelta = cellfun(derv, atTimes(cfp,:), fnOpts{:}, erOpts{:});
        itDelta = cellfun(derv, itTimes(:,cfp), fnOpts{:}, erOpts{:});
        ddm = cellfun(dsmt,...
            itDelta(:, cfp), flip(atDelta(cfp,:)'), fnOpts{:});
        dm = cellfun(@(x,y) y - x(:,1)', itTimes(:,cfp), flip(atTimes(cfp,:)'),...
            fnOpts{:});
    end
% Function to eliminate a trigger that was detected in the roller position
% file but didn't actually happen according to the ground truth: the intan
% trigger signals
    function detectFalseAlarms()
        itTimes(:,cfp) = cellfun(@(x) x./3e4, tSubs, fnOpts{:});
        % Registering and deleting the empty triggers
        eiFlag = cellfun(@(x) isempty(x), itTimes); itTimes(eiFlag) = [];
        if size(itTimes,1) ~= size(atTimes,2)
            atTimes(flip(eiFlag)) = [];
            fprintf(1, "%s was a false detection!\n", atNames(flip(eiFlag)))
            atNames(flip(eiFlag)) = [];
        end
    end
% Read the trigger file and compute the triggers directly without keeping
% the trigger signals in memory
    function subs = readTriggerFile(tfName)
        fID = fopen(tfName,'r'); trig = fread(fID, [2,Inf], 'uint16=>int32');
        [~] = fclose(fID); trig = int16(trig - median(trig, 2));
        trig(1,:) = -trig(1,:);
        subs = getSubsFromTriggers(trig);
    end
% Detection of stimulation pairs between arduino and intan
    function correctAnomalities()
        [Nta, Nti] = cellfun(@size, dm);
        errTh = -3; 
        for cc = 1:size(dm)
            ca = cc;
            if flipFlag
                ca = 3-cc;
            end
            % Taking a guess from the number of triggers in both trigger
            % recorders
            if Nti(cc) > Nta(cc)
                fprintf(1, "Perhaps the arduino missed N trigger(s)\n");
            elseif Nta(cc) > Nti(cc)
                fprintf(1, "Perhaps the arduino added N trigger(s) from noise\n");
            else
                fprintf(1, "Seems like the triggers have good pairing\n");
                intSim = diag(log2(abs(ddm{cc}+1))); mm = mean(intSim);
                sm = std(intSim); ffm = sm/mm;
                if abs(ffm) > 0.8
                    % Perhaps there's a change in delay or one trigger is
                    % slightly more shifted than the others
                    fprintf(1, "Viewing the triggers recommended\n")
                    fprintf(1, "Perhaps need to reposition a trigger\n")
                else
                    % Seems like there's no problem at all ;)
                    fprintf(1, "Looks like this set of triggers have no issue\n")
                    continue
                end
            end
            preArd = dm{cc} < 0; [~, mnSubs] = sort(dm{cc}(:), "ascend",...
                "ComparisonMethod", "abs");
            mnSubs(~preArd(mnSubs)) = [];
            [aSubs, iSubs] = ind2sub(size(dm{cc}), mnSubs);
            mxSub = min(Nti(cc),Nta(cc));
            cp = 1; r2 = 0; err = Inf;
            while err > errTh || r2 < 0.98
                while length(unique(iSubs(cp:mxSub-1+cp))) < mxSub ||...
                        length(unique(aSubs(cp:mxSub-1+cp))) < mxSub
                    cp = cp + 1;
                end
                xSubs = sort(iSubs(cp:mxSub-1+cp));
                x = itTimes{cc}(xSubs,1);
                ySubs = sort(aSubs(cp:mxSub-1+cp));
                y = atTimes{ca}(ySubs);
                [mdl, yhat, r2] = fit_poly(x, y, 1);
                err = log(norm(y - yhat, 1));
                cp = cp + 1;
            end
            if Nta(cc) > Nti(cc)
                % More arduino pulses
                atTimes{ca}(setdiff(1:Nta(cc), ySubs)) = [];
            elseif Nta(cc) < Nti(cc)
                % Missed arduino pulses
                atTimesHat = (itTimes{cc}(setdiff(1:Nti(cc), xSubs),1).^[1,0])*mdl;
                atTimes2 = zeros(size(itTimes{cc},1),1);
                atTimes2(ySubs) = atTimes{ca}(:);
                atTimes2(setdiff(1:Nti(cc), ySubs)) = atTimesHat;
                atTimes{ca} = atTimes2;
            else
                % Correct cardinality
            end
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
Nrp = numel(rpFiles); Nit = numel(itFiles);
if Nrp ~= Nit
    fprintf(1, "Discrepancy between files!\n")
    fprintf(1, "# roller files: %d, # trigger files: %d\n",...
        numel(rpFiles), numel(itFiles));
    return
end
% Extract the date from the file names. Do the names correspond?
rpDates = getDates(rpFiles, "Roller_position");
itDates = getDates(itFiles, "TriggerSignals");
if any((rpDates - itDates) ~= 0)
    fprintf(1, "The file names do not correspond!\n")
    fprintf(1, "Please check and maybe even organise the files in this ")
    fprintf(1, "folder\n"); fprintf(1, "Roller files:\n");
    fprintf(1, "%s\n", arrayfun(@(x) string(x.name), rpFiles))
    fprintf(1, "Trigger files:\n");
    fprintf(1, "%s\n", arrayfun(@(x) string(x.name), itFiles))
    return
end

%% Searching for arduino instabilities
itNames = ["P", "L"];

for cfp = 1:numel(rpFiles)
    flipFlag = false;
    % Read file for arduino trigger times and names
    [atTimes, atNames] = getArduinoTriggers(...
        fullfile(rpFiles(cfp).folder, rpFiles(cfp).name));
    if ~all(itNames == atNames)
        flipFlag = true;
    end
    % Read trigger signals
    tfName = fullfile(itFiles(cfp).folder,itFiles(cfp).name);
    % Subs from the trigger signals
    tSubs = readTriggerFile(tfName); detectFalseAlarms()
    % Computing the similarity between the trigger intervals from arduino
    % and intan
    [ddm, dm] = computeSimilarityMatrices();
    correctAnomalities();
end
end

function A = falseLaserDetection(S, varargin)
fprintf(1, "No laser detection!\n")
fprintf(1, "%s\n", S.message)
A = [];
end
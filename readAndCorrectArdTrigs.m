function [atTimes, atNames, itTimes, itNames, Nt] = readAndCorrectArdTrigs(dataDir)
%READANDCORRECTARDTRIGS pairs the trigger times from the arduino board with
%intan onset trigger times and saves the times in seconds in a file called
%ArduinoTriggersYYYY-MM-DDTHH_mm_ss.mat accordingly with the
%Roller_position and TriggerSignal files with the same date in the file
%name.
%   Detailed explanation goes here
%% Auxiliary variables and functions
fnOpts = {"UniformOutput", false};
dateFormStr = 'yyyy-MM-dd''T''HH_mm_ss';
outFileFormat = "ArduinoTriggers%s.mat";
derv = @(x) diff(x(:,1));
dsmt = @(x,y) distmatrix(x, y, 2);
erOpts = {"ErrorHandler", @falseLaserDetection};

fs = 3e4;
% Function to get the edges from the trigger signals in the trigger file.
    function tSubs = getSubsFromTriggers(trig)
        % Unstable line! Display function not working!
        tObj = arrayfun(@(x) StepWaveform(trig(x,:), 3e4, ...
            'verbose',false), [1;2]);
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
        atTimes = cell(1,1); atNames = "";
        % Read the roller position file
        [rp, tTimes] = readRollerPositionsFile(rpFilePath);
        % Converting the trigger times from microseconds to seconds
        if ~isempty(tTimes)
            tTimes(1,:) = cellfun(@(x) (x - rp(1,2))/1e6, tTimes(1,:), fnOpts{:});
            % Deleting noise entries to the arduino (times with less than 1 second
            % difference)
            delFlags = cellfun(@(x) [false;diff(x) < 1], tTimes(1,:), fnOpts{:});
            tTimes(1,:) = cellfun(@(x,y) x(~y), tTimes(1,:), delFlags, fnOpts{:});
            % Arduino times
            atTimes = tTimes(1,:); atNames = string(tTimes(2,:));
        end
    end
% Function to eliminate a trigger that was detected in the roller position
% file but didn't actually happen according to the ground truth: the intan
% trigger signals
    function itTimes = detectFalseAlarms()
        itTimes = cellfun(@(x) x./3e4, tSubs, fnOpts{:});
        % Registering and deleting the empty triggers
        eiFlag = cellfun(@(x) isempty(x), itTimes); itTimes(eiFlag) = [];
        if size(itTimes,1) ~= size(atTimes,2)
            atTimes(flip(eiFlag)) = [];
            fprintf(1, "%s was a false detection!\n", atNames(flip(eiFlag)))
            atNames(flip(eiFlag)) = [];
            itNames(eiFlag) = [];
        end
    end
% Function to compute 2 matrices: similarity between trigger intervals
% (ddm) and similarity between trigger times (dm)
    function [ddm, dm] = computeSimilarityMatrices()
        atDelta = cellfun(derv, atTimes, fnOpts{:}, erOpts{:});
        itDelta = cellfun(derv, itTimes, fnOpts{:}, erOpts{:});
        ddm = cellfun(dsmt, itDelta, flip(atDelta'), fnOpts{:});
        dm = cellfun(@(x,y) y-x(:,1)', itTimes, flip(atTimes'),fnOpts{:});
    end
% Read the trigger file and compute the triggers directly without keeping
% the trigger signals in memory
    function subs = readTriggerFile(tfName)
        fID = fopen(tfName,'r'); trig = fread(fID, [2,Inf], 'uint16=>int32');
        [~] = fclose(fID); trig = int16(trig - median(trig, 2));
        trig(1,:) = -trig(1,:); Ns = length(trig);
        subs = getSubsFromTriggers(trig);
    end
% Detection of stimulation pairs between arduino and intan
    function correctAnomalities()
        [Nta, Nti] = cellfun(@size, dm);
        % errTh = -3;
        Ntrigs = length(dm);
        for cc = 1:Ntrigs
            ca = cc;
            if flipFlag && Ntrigs > 1
                ca = 3-cc;
            end
            % Taking a guess from the number of triggers in both trigger
            % recorders
            if Nti(cc) > Nta(cc)
                fprintf(1, "Perhaps the arduino missed N trigger(s)\n");
            elseif Nta(cc) > Nti(cc)
                fprintf(1, "Perhaps the arduino added N trigger(s) from noise\n");
            else
                intSim = diag(log2(abs(ddm{cc}+1))); mm = mean(intSim);
                sm = std(intSim); ffm = sm/mm;
                if abs(ffm) > 1
                    % Perhaps there's a change in delay or one trigger is
                    % slightly more shifted than the others
                    fprintf(1, "Viewing the triggers recommended\n")
                    fprintf(1, "Perhaps need to reposition a trigger\n")
                    fprintf(1, "Or it is just a change in delay\n")
                    continue
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
            mxSub = min(Nti(cc),Nta(cc)); naiveSubs = 1:Nti(cc);
            pairFlags = abs(iSubs - aSubs) < abs(Nti(cc) - Nta(cc)) + 1;
            dstPrs = dm{cc}(mnSubs(pairFlags)); 
            unqePrs = unique([iSubs(pairFlags), aSubs(pairFlags), dstPrs], ...
                "row");
            % Correlation between trigger times
            [sbXc, aLag] = xcorr(atTimes{ca}, itTimes{cc}(:,1));
            [~, mxLag] = max(sbXc); aLag = aLag(mxLag);
            dirFlag = diff(unqePrs(:,1:2), 1, 2) == aLag:-sign(aLag):0;
            dirFlag = any(dirFlag, 2);  %#ok<NASGU> 
            % Removing pairs going in the opposite direcion of the
            % cross-corrlation maximum lag. ???
            % % % % % unqePrs(~dirFlag,:) = [];
            % Distribution of the distances found between the likely pairs
            dstDom = linspace(min(unqePrs(:,3))*1.1, 0.1, 257)';
            % Quartiles
            dstDist = ksdensity(dstPrs, dstDom, "Bandwidth", ...
                1/sqrt(2));
            dstCnt = mean([dstDom(1:end-1),dstDom(2:end)], 2);
%             Q = arrayfun(@(q) fminbnd(@(x) ...
%                 norm(interp1(dstDom, cumsum(dstDist)-q, x),2), ...
%                 dstDom(1), dstDom(end)), 1/4 * [1,3]);
            Q = quantile(dstPrs, 1/4 * [1,3]);
            modDst = Q(1)*1.1;
            % Log-dist
            % [bC, bE] = prepareLogBinEdges(dstPrs, 64);
            % bN = histcounts(log10(abs(dstPrs)), bE, 'Normalization','pdf');
            % Kernel distribution
            dstDist2 = ksdensity(dstPrs, dstDom, "Bandwidth", ...
                floor(range(dstPrs)/sqrt(sum(pairFlags)*pi)));
            % Estimation for the most repeated delay (the actual delay
            % between intan and arduino)
            loPks = arrayfun(@(x) fminsearch(@(y) -interp1(dstDom, ...
                dstDist2, y), x), rand(100,1)*(min(dstPrs)*1.1));
            dstPks = [uniquetol(loPks, 0.01/max(abs(loPks)));...
                -mode(abs(dstPrs))];
            % Likelihood of the pair given the distance distribution
%             lk = interp1(dstDom, dstDist2*(1./max(dstDist2)), ...
%                 unqePrs(:,3), "pchip", "extrap");
            lk = interp1(dstDom, dstDist*(1./max(dstDist)), unqePrs(:,3), ...
                "pchip", "extrap");
            xyM = zeros(Nti(cc),3);
            for ci = naiveSubs
                iFlag = unqePrs(:,1) == ci;
                if sum(iFlag)
                    [~, iSb] = max(lk(iFlag));
                    xyM(ci,:) = unqePrs(find(cumsum(iFlag) == iSb, 1, ...
                        "first"),:);
                end
            end
            xyM(~any(xyM,2),:) = [];
            unqePrs = xyM;
            % Initialization of the logical reduction loop
            allFlag = unqePrs(:,1:2) == reshape(naiveSubs, 1,1,[]);
            IoA = all(any(allFlag,3))'; snglPairs = zeros(mxSub, 2);
            if all(IoA)
                [~, IoASub] = min([Nti(cc),Nta(cc)]); IoA = false(2,1);
                IoA(IoASub) = true;
            end
            scndFlag = false;
            auxCnt = 1;
            % Logical reduction of possibilities
            while auxCnt <= mxSub && ~isempty(unqePrs)
                allFlag = unqePrs(:,1:2) == reshape(naiveSubs, 1,1,[]);
                appCount = squeeze(sum(allFlag, 1));
                appCount(appCount == 0) = NaN;
                % Lower cardinality trigger times with only one pairing possibility
                mstUseFlag = any(appCount == IoA);
                if all(~mstUseFlag)
                    unlikelyFlag = isnan(appCount(IoA,:));
                    unpairFlag = any(unqePrs(:,IoA) == naiveSubs(unlikelyFlag),2);
                    if all(~unpairFlag)
                        unlikelyFlag = unqePrs(:,3) < modDst*1.1;
                        if all(~unlikelyFlag)
                            fprintf(1, "Taking the pair(s) that is(are) closest to")
                            fprintf(1, " the linear fit\n")
                            xSubs = snglPairs(snglPairs(:,1) ~= 0,1);
                            ySubs = snglPairs(snglPairs(:,2) ~= 0,2);
                            x = itTimes{cc}(xSubs,1);
                            y = atTimes{ca}(ySubs);
                            [n, d] = getHesseLineForm(fit_poly(x, y, 1));
                            M = [itTimes{cc}(unqePrs(:,1),1),...
                                atTimes{ca}(unqePrs(:,2))];
                            yErr = log(abs(M*n - d));
                            bigErrFlag = yErr > 0;
                            if all(~bigErrFlag) || all(bigErrFlag)
                                [~, furthSub] = max(yErr);
                                unqePrs(furthSub,:) = [];
                            else
                                unqePrs(bigErrFlag,:) = [];
                            end
                        else
                            fprintf(1, "Removing pairs with a 'big'")
                            fprintf(1, " distance. (%.2f)\n", modDst*1.1)
                            unqePrs(unlikelyFlag,:) = [];
                        end
                    else
                        unqePrs(unpairFlag,:) = [];
                    end
                    continue
                end
                % Pairs that *must* be in the pair selection
                unPrFlag = any(unqePrs(:,IoA) == naiveSubs(mstUseFlag),2);
                Nup = sum(unPrFlag);
                snglPairs(auxCnt:auxCnt-1+Nup,:) = unqePrs(unPrFlag,1:2);
                auxCnt = auxCnt+Nup;
                unavFlags = arrayfun(@(x) any(unqePrs(:,x) == snglPairs(:,x)',2), 1:2,...
                    "UniformOutput", false); unavFlags = [unavFlags{:}];
                unqePrs(any(unavFlags,2),:) = [];
                if isempty(unqePrs) && ~scndFlag
                    if nnz(snglPairs) < numel(snglPairs)
                        snglPairs = zeros(mxSub, 2); auxCnt = 1;
                        fprintf(1, "There was a miscalculation in the ")
                        fprintf(1, "trigger assignment!\n")
                        unqePrs = unique([iSubs(pairFlags), ...
                            aSubs(pairFlags), dstPrs], "row");
                        ddstPrs = unqePrs(:,3) - dstPks';
                        prFlag = any(ddstPrs < 1 & ddstPrs > -1,2);
                        allFlag = unqePrs(prFlag,1:2) == reshape(naiveSubs, 1,1,[]);
                        appCount = squeeze(sum(allFlag, 1));
                        appCount(appCount == 0) = nan;
                        mstUseFlag = any(appCount == IoA);
                        auxUP = unqePrs(prFlag, :);
                        unPrFlag = any(unqePrs(prFlag,IoA) == naiveSubs(mstUseFlag),2);
                        Nup = sum(unPrFlag);
                        snglPairs(auxCnt:auxCnt-1+Nup,:) = auxUP(unPrFlag,1:2);
                        auxCnt = auxCnt+Nup;
                        unavFlags = arrayfun(@(x) any(unqePrs(:,x) == snglPairs(:,x)',2), 1:2,...
                            "UniformOutput", false); unavFlags = [unavFlags{:}];
                        unqePrs(any(unavFlags,2),:) = [];
                        scndFlag = true;
                        continue
                    end
                end
            end
            snglPairs(~any(snglPairs,2),:) = [];
            snglPairs = sortrows(snglPairs, 1);
            raFlag = [1;diff(snglPairs(:,2))]<1; % repeated arduino trigger
            snglPairs(raFlag,:) = [];
            x = itTimes{cc}(snglPairs(:,1),1);
            y = atTimes{ca}(snglPairs(:,2));
            [mdl, ~, r2] = fit_poly(x, y, 1);
            % [mdl, yhat, r2] = fit_poly(x, y, 1);
            % err = log(norm(y - yhat,1));
            if r2 < 0.8
                fprintf(1, "Better take a look at the triggers!\n")
                return
            end
            mssSubs = setdiff(1:Nti(cc), snglPairs(:,1));
            M = (itTimes{cc}(mssSubs,1).^[1,0]);
            atTimesHat = M*mdl;
            atTimes2 = zeros(Nti(cc),1);
            atTimes2(snglPairs(:,1)) = atTimes{ca}(snglPairs(:,2));
            atTimes2(mssSubs) = atTimesHat;
            atTimes{ca} = atTimes2;
            %{
            if Nta(cc) > Nti(cc)
                % More arduino pulses
                atTimes{ca}(setdiff(1:Nta(cc), snglPairs(:,2))) = [];
            elseif Nta(cc) < Nti(cc)
                % Missed arduino pulses
                mssSubs = setdiff(1:Nti(cc), snglPairs(:,1));
                M = (itTimes{cc}(mssSubs,1).^[1,0]);
                %{ 
                % We assume that we miss only 1 trigger. 
                if err > errTh
                    xSubs = snglPairs(snglPairs(:,1)>mssSubs(1),1);
                    ySubs = snglPairs(snglPairs(:,1)>mssSubs(1),2);
                    x = itTimes{cc}(xSubs,1);
                    y = atTimes{ca}(ySubs);
                    mdl = fit_poly(x, y, 1);
                end
                %}
                atTimesHat = M*mdl;
                atTimes2 = zeros(Nti(cc),1);
                atTimes2(snglPairs(:,1)) = atTimes{ca};
                atTimes2(mssSubs) = atTimesHat;
                atTimes{ca} = atTimes2;
            else
                % Correct cardinality
            end
            %}
        end
    end
%% Validation section
% Given folder exists?
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
Ns = 0;
%% Searching for arduino instabilities
var2save = {'atTimes', 'atNames', 'itTimes', 'itNames', 'Nt', 'minOfSt'};
var2Check = var2save;
flfl = @(x, y) fullfile(x, y);
fobPthIdx = @(x, y) fullfile(x(y).folder, x(y).name);
fobPth = @(x) fullfile(x.folder, x.name);
lcFiles = dir(flfl(string(dataDir), "Laser*.csv"));
pcFiles = dir(flfl(string(dataDir), "Puff*.csv"));
for cfp = 1:numel(rpFiles)
    itNames = ["P"; "L"];
    flipFlag = false;
    atFileName = flfl(rpFiles(cfp).folder,...
        sprintf(outFileFormat, string(rpDates(cfp), dateFormStr)));
    atMObj = matfile(atFileName); varIn = who(atMObj);
    if ~exist(atFileName, "file") || ~all(contains(varIn,var2Check))
        % Read file for *arduino* trigger times and names
        [atTimes, atNames] = getArduinoTriggers(fobPthIdx(rpFiles, cfp));
        % Read trigger signals and extracet subs from the *intan* trigger
        % signals
        tfName = fobPthIdx(itFiles, cfp); tSubs = readTriggerFile(tfName);
        if isempty(atTimes{1})
            Nti = cellfun(@(x) size(x, 1), tSubs(:));
            % Old version without trigger times in the roller position
            trigFls = [lcFiles(cfp), pcFiles(cfp)];
            [atTimes, atNames] = arrayfun(@getCSVTriggers,...
                arrayfun(@(x) string(fobPth(x)),trigFls), fnOpts{:});
            atNames = cat(2, atNames{:});
            [iS, aS, itNames] = matchOrder(atNames, itNames);
            Nta = cellfun(@(x) size(x, 1), atTimes(:)); 
            trm = Nta(aS) - Nti(iS);
            atTimes = cellfun(@(x,y) x(y+1:end,1) - (10 + ...
                second(rpDates(cfp), "secondofday")), atTimes, ...
                num2cell(trm)', fnOpts{:});
            itTimes = detectFalseAlarms(); Nt = Ns/fs;
        else
            if ~all(itNames == atNames)
                flipFlag = true;
            end
            itTimes = detectFalseAlarms(); Nt = Ns/fs;
            % Computing the similarity between the trigger intervals from
            % arduino and intan
            [ddm, dm] = computeSimilarityMatrices();
            correctAnomalities();
            [iS, aS, itNames] = matchOrder(atNames, itNames); 
        end
        ofSts = arrayfun(@(x) itTimes{iS(x)}(1,1) - atTimes{aS(x)}(1), ...
            1:size(atNames)); minOfSt = min(ofSts);
        fprintf(1, "Correcting arduino triggers by %.3f seconds\n", minOfSt)
        atTimes = cellfun(@(x) x+minOfSt, atTimes, fnOpts{:});
        fprintf(1, "Saving %s...\n", atFileName)
        save(atFileName, var2save{:});
    else
        fprintf(1, "File already saved! Skipping...\n")
    end
end

end

function [iS, aS, itNames] = matchOrder(atNames, itNames)
tmMat = atNames == itNames;
itNames(~any(tmMat,2)) = []; tmMat(~any(tmMat,2)) = [];
[iS, aS] = find(tmMat);
end

function A = falseLaserDetection(S, varargin)
fprintf(1, "No laser detection!\n")
fprintf(1, "%s\n", S.message)
A = [];
end
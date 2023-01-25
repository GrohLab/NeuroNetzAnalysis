function [resStruct] = analyseBehaviour(behDir, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Auxiliary variables
% Input files
dlcPttrn = 'roller*filtered.csv';
% rpPttrn = "Roller_position*.csv"; vdPttrn = "roller*.avi";
% Output files
afPttrn = "ArduinoTriggers*.mat"; rfPttrn = "RollerSpeed*.mat";
% DLC variables
dlcVars = {'behDLCSignals', 'dlcNames'};
% lPttrn = "Laser*.csv"; pPttrn = "Puff*.csv";
% Functional cell array
fnOpts = {'UniformOutput', false}; 
% Aesthetic cell arrays
axOpts = {'Box','off','Color','none', 'NextPlot', 'add'};
lgOpts = cat(2, axOpts{1:2}, {'Location','best'});
ptOpts = {"Color", 0.7*ones(1,3), "LineWidth", 0.2;...
    "Color", "k", "LineWidth",  1.5};
phOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.25, 'FaceColor'};
% Anonymus functions
behHere = @(x) fullfile(behDir, x);
flfile = @(x) fullfile(x.folder, x.name); iemty = @(x) isempty(x);
istxt = @(x) isstring(x) | ischar(x);
mat2ptch = @(x) [x(1:end,:)*[1;1]; x(end:-1:1,:)*[1;-1]];
getSEM = @(x, idx) [mean(x(:,idx),2), std(x(:,idx),1,2)./sqrt(sum(idx))];
% getQs = @(x, idx) quantile(x(:,idx),3,2);
% q2patch = @(x) [x(:,1);x(end:-1:1,3)];
getThreshCross = @(x) sum(x)/size(x,1);
m = 1e-3; k = 1/m;
%% Input validation
p = inputParser;

defVW = [-250, 500]*m; %milliseconds
defRW = [5, 400]*m;
checkVW = @(x) all([isnumeric(x), numel(x) == 2, x(1)<x(2)]);

% Puff or laser (or whatever these would represent)
defCond = "P"; % Default condition to center the behavioural signals
checkCond = @(x) numel(x) == 1 & istxt(x) & (strcmpi(x,"P") | strcmpi(x,"L"));

% If the triggers don't represent the same for all the experiment or if you
% want to see only a subset of all the triggers.
defTrigSubs = "all";
checkTrigSubs = @(x) strcmpi(x,"all") | isPositiveIntegerValuedNumeric(x);
% Check paired stimulus
checkPF = @(x) islogical(x) & (ismatrix(x) | isvector(x));
% Check conditions names
checkCondNames = @(x) (iscell(x) | istxt(x)) & numel(x) >= 1;

defSpeedTh = 0.1:0.1:3; % Speed thresholds
checkSpeedTh = @(x) all([diff(x) > 0, x > 0, numel(x) > 1]);

checkTH = @(x) x > 0 & numel(x) == 1;
defSig = 2.5; % Spontaneous standard deviation
defSMed = 0.2; % Spontaneous Median
defMed = 1; % Median in all viewing window

addRequired(p, 'behDir', istxt)
addParameter(p, 'ViewingWindow', defVW, checkVW)
addParameter(p, 'ResponseWindow', defRW, checkVW)
addParameter(p, 'Condition', defCond, checkCond)
addParameter(p, 'TriggerSubscripts', defTrigSubs, checkTrigSubs)
addParameter(p, 'PairedFlags', 'none', checkPF)
addParameter(p, 'ConditionsNames', 'none', checkCondNames)
addParameter(p, 'SpeedTh', defSpeedTh, checkSpeedTh)
addParameter(p, 'SponSigTh', defSig, checkTH)
addParameter(p, 'SponMedianTh', defSMed, checkTH)
addParameter(p, 'MedianTh', defMed, checkTH)
addParameter(p, 'FigureDirectory', "Figures", istxt)
addParameter(p, 'verbose', true, @(x) islogical(x) & numel(x) == 1)

parse(p, behDir, varargin{:})

behDir = p.Results.behDir;
bvWin = p.Results.ViewingWindow;
brWin = p.Results.ResponseWindow;
chCond = p.Results.Condition;
trigSubs = p.Results.TriggerSubscripts;
pairedStim = p.Results.PairedFlags;
consCondNames = {p.Results.ConditionsNames};
spTh = {p.Results.SpeedTh};
sigTh = p.Results.SponSigTh;
sMedTh = p.Results.SponMedianTh;
tMedTh = p.Results.MedianTh;
figureDir = p.Results.FigureDirectory;
verbose = p.Results.verbose;

if bvWin(1) > brWin(1) || bvWin(2) < brWin(2)
    if verbose
        fprintf(1, 'Given response window lays outside the viewing window!\n')
        fprintf(1, 'VW:%.2f - %.2f ms\nRW:%.2f - %.2f ms\n', bvWin*k, brWin*k)
        fprintf(1, 'Resetting to default windows!\n')
    end
    bvWin = defVW; brWin = defRW;
    if verbose
        fprintf(1, 'VW:%.2f - %.2f ms\nRW:%.2f - %.2f ms\n', bvWin*k, brWin*k)
    end
end

if ~exist(behDir, "dir")
    if verbose
        fprintf(1, 'Given folder doesn''t exist!\nExiting...')
    end
    return
end

awConfStruct = struct('Condition', chCond, 'Triggers', trigSubs, ...
    'ViewingWindow', bvWin, 'ResponseWindow', brWin, 'SpeedTh', spTh, ...
    'Thresholds', struct('SpontSigma', sigTh, 'SpontMedian', sMedTh, ...
    'MedianTh', tMedTh)); %#ok<NASGU> 

% Validation of the figure directory
[figParentDir, ~] = fileparts(figureDir);
if ~strlength(figParentDir)
    figureDir = behHere(figureDir);
    if ~exist(figureDir, "dir")
        if verbose
            fprintf(1, "Creating figure directory.\n")
        end
        if ~mkdir(figureDir)
            if verbose
                fprintf(1, "Unable to create figure directory!\n")
            end
        end
    end
elseif ~exist(figureDir, "dir")
    fprintf(1, "Given figure directory doesn't exist!")
    fprintf(1, "Creating 'Figures' in the given directory.\n")
    figureDir = behHere("Figures");
    if ~exist(figureDir, "dir")
        if ~mkdir(figureDir) && verbose
            fprintf(1, "Unable to create figure directory!\n")
        end
    end
end

%% Create or read output files
if verbose
    fprintf(1,'Time window: %.2f - %.2f ms\n',bvWin*k)
    fprintf(1,'Response window: %.2f - %.2f ms\n',brWin*k)
end

if iemty(dir(behHere(afPttrn)))
    readAndCorrectArdTrigs(behDir);
end
% Roller speed
rfFiles = dir(behHere(rfPttrn));
if iemty(rfFiles)
    [~, vf, ~, fr, Texp] = createRollerSpeed(behDir);
    rfFiles = dir(behHere(rfPttrn));
end
if numel(rfFiles) == 1
    rfName = fullfile(rfFiles.folder, rfFiles.name);
    load(rfName) %#ok<LOAD>
    try
        %              Encoder steps  Radius^2
        en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*fr;
    catch
        try
            %              Encoder steps  Radius^2
            en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*rollFs;
            fr = rollFs;
        catch
            en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*fsRoll;
            fr = fsRoll;
        end
    end
end

% DLC signals
dlcFiles = dir(behHere(dlcPttrn));
if ~iemty(dlcFiles)
    behPttrn = "BehaveSignals%s%s.mat";
    endng = "DLC" + string(extractBetween(dlcFiles(1).name, "DLC", ".csv"));
    behPath = joinBehDates(dlcFiles, "roller", behHere(behPttrn), ...
        endng);
    if ~exist(behPath, "file")
        if verbose
            fprintf(1, "Computing signals... \n")
        end
        dlcTables = arrayfun(@(x) readDLCData(flfile(x)), dlcFiles, fnOpts{:});
        a_bodyParts = cellfun(@(x) getBehaviourSignals(x), ...
            dlcTables, fnOpts{:});
        behStruct = cellfun(@(x) getWhiskANoseFromTable(x), a_bodyParts);
        % Butter order 3 low-pass filter 35 Hz cut off frequency @ 3 dB
        [b, a] = butter(3, 70/fr, 'low');
        behDLCSignals = cat(1, behStruct(:).Signals);
        dlcNames = behStruct(1).Names; clearvars behStruct;
        behDLCSignals = filtfilt(b, a, behDLCSignals);
        save(behPath, dlcVars{:})
        if verbose
            fprintf(1, "Saving... Done!\n")
        end
    else
        if verbose
            [~, auxFile] = fileparts(behPath);
            fprintf(1, "Loading %s\n", auxFile)
        end
        dlcFileVars = load(behPath, dlcVars{:});
        behDLCSignals = dlcFileVars.behDLCSignals;
        dlcNames = dlcFileVars.dlcNames; clearvars dlcFileVars;
    end
else
    if verbose
        fprintf(1, "No filtered .csv DLC files found!\n")
    end
    behDLCSignals = [];
    dlcNames = "";
end
Nbs = size(behDLCSignals, 2);

% Triggers
atVar = {'atTimes', 'atNames', 'itTimes', 'itNames'};
afFiles = dir(behHere(afPttrn));
if ~iemty(afFiles)
    jointArdPttrn = "JointArduinoTriggers%s%s.mat";
    jointArdOut = behHere(jointArdPttrn);
    jointArdPath = joinBehDates(afFiles, extractBefore(afPttrn,"*"), ...
        jointArdOut);
    [~, auxFile] = fileparts(jointArdPath);
    if exist(jointArdPath, "file")
        if verbose
            fprintf(1, "Loading %s... ", auxFile)
        end
        atV = load(jointArdPath, atVar{:});
        atTimes = atV.atTimes;
        atNames = atV.atNames;
        itTimes = atV.itTimes;
        itNames = atV.itNames; clearvars atV;
        if verbose
            fprintf(1, "done!\n")
        end
    else
        % Arduino
        atV = arrayfun(@(x) load(flfile(x), atVar{:}), afFiles);
        Nrecs = length(atV);
        atT = arrayfun(@(x, z) cellfun(@(y, a) y+a, ...
            x.atTimes, repmat(z,1,length(x.atTimes)), fnOpts{:}), atV', ...
            num2cell(cumsum([0, Texp(1:end-1)])), fnOpts{:});
        trig_per_recording = cellfun(@(x) size(x,2), atT);
        [max_trigs, record_most_trigs] = max(trig_per_recording);
        record_trig_cont_ID = arrayfun(@(x) ...
            contains(atV(record_most_trigs).atNames, x.atNames), ...
            atV(1:Nrecs), fnOpts{:}); outCell = cell(Nrecs, max_trigs);
        for cr = 1:Nrecs
            outCell(cr,record_trig_cont_ID{cr}) = atT{cr};
        end
        atTimes = arrayfun(@(x) cat(1, outCell{:,x}), 1:size(outCell,2), ...
            fnOpts{:});
        atNames = atV(1).atNames;
        %{
        atT = cat(1, atT{:});
        atTimes = arrayfun(@(x) cat(1, atT{:,x}), 1:size(atT,2), fnOpts{:});
        atNames = atV(1).atNames;
        %}
        % Intan
        itT = arrayfun(@(x, z) cellfun(@(y, a) y+a, ...
            x.itTimes, repmat(z,size(x.itTimes)), fnOpts{:}), atV', ...
            num2cell(cumsum([0, Texp(1:end-1)])), fnOpts{:});
        outCell = cell(Nrecs, max_trigs);
        for cr = 1:Nrecs
            outCell(cr,record_trig_cont_ID{cr}) = itT{cr};
        end
        itTimes = arrayfun(@(x) cat(1, outCell{:,x}), 1:size(outCell,2), ...
            fnOpts{:});
        itNames = atV(1).itNames;
        if numel(afFiles) > 1
            if verbose
                fprintf(1, "Saving %s...\n", auxFile)
            end
            save(jointArdPath, atVar{:})
        end
    end
end

%% Cutting the behavioural signals
aSub = strcmpi(atNames, chCond);
iSub = strcmpi(itNames, chCond);
if isnumeric(trigSubs)
    Ntrig = size(atTimes{aSub},1);
    triggPassFlag = Ntrig < trigSubs;
    if any(triggPassFlag)
        if verbose
            fprintf(1, "Given subscripts pass the maximum trigger number")
            fprintf(1, " %d!\n", Ntrig)
            fprintf(1, "Clipping the given subscripts to the available ")
            fprintf(1, "range!\n")
            fprintf(1, "Deleting: %d", trigSubs(triggPassFlag))
            fprintf(1, "/n")
        end
        trigSubs(triggPassFlag) = [];
    end
elseif istxt(trigSubs) && strcmpi(trigSubs, "all")
    trigSubs = 1:size(atTimes{aSub},1);
end
% Roller speed
[~, vStack] = getStacks(false, round(atTimes{aSub}(trigSubs)*fr), 'on', bvWin,...
    fr, fr, [], vf*en2cm); [Nbt, Ntr] = size(vStack, [2, 3]); Na = Ntr;
% Connection to DE_Jittering or as a stand-alone function.
if istxt(pairedStim) && strcmpi(pairedStim, "none")
    % Stand alone function without stimulus pairing
    pairedStim = true(Na, 1);
    Nccond = 1;
    consCondNames = cellstr(atNames(aSub));
elseif islogical(pairedStim)
    % Connected to DE_Jittering
    NTa = sum(pairedStim(:));
    if size(pairedStim, 1) ~= Ntr
        fprintf(1, "The number of trials in the given flags (%d) do not", ...
            NTr)
        fprintf(1, " correspond to the number of trials in behaviour files");
        fprintf(1, " (%d)!\n", Ntr);
    else
        Nccond = size(pairedStim, 2);
        if size(pairedStim,1)~=Ntr
            pairedStim(all(~pairedStim,2),:) = [];
        end
        while numel(consCondNames) ~= Nccond || all(~cellfun(istxt, consCondNames))
            consCondNames = consCondNames{:};
            if ~iscell(consCondNames)
                fprintf(1, "Considered condition names causing trouble!\n")
                disp(consCondNames)
                fprintf(1, "Naming the conditions with a consecutive letter\n")
                % 'a' has ASCII code 97
                consCondNames = arrayfun(@(x) sprintf('%s', x), ...
                    97:97+Nccond-1, fnOpts{:});
                disp(consCondNames)
                break
            end
        end
    end
end
% Estimating the setpoint for whiskers and nose signals.
behTx_all = (0:size(behDLCSignals, 1) - 1)/fr;
spFile = fullfile(behDir, "SetPoint.mat");
if ~exist(spFile, "file")
    set_point = arrayfun(@(x) fitSpline(behTx_all, behDLCSignals(:,x), 1, ...
        0.3, 0.95, verbose), 1:size(behDLCSignals,2), fnOpts{:});
    set_point = [set_point{:}];
    if ~isempty(set_point)
        save(spFile, "set_point")
    end
else
    set_point = load(spFile,"set_point"); set_point = set_point.set_point;
end

% Whiskers and nose
[~, dlcStack] = getStacks(false, round(itTimes{iSub}(trigSubs,:)*fr), 'on', ...
    bvWin, fr, fr, [], [behDLCSignals, set_point]');
behStack = cat(1, dlcStack(1:Nbs,:,:), vStack);
Nbs = size(behStack, 1);
behStack = arrayfun(@(x) squeeze(behStack(x, :, :)), (1:size(behStack,1))',...
    fnOpts{:});
behNames = [dlcNames, "Roller speed"]; behNames(~strlength(behNames))=[];
rollYL = "Roller speed [cm/s]";
if size(behNames,2)>1
    yLabels = [repmat("Angle [°]", 1, 3), rollYL];
else
    yLabels = rollYL;
end

% Whiskers, nose, and roller speed

%TODO: include user input for the thresholds for whiskers and nose signals.
if Nbs>1
    % Probability thresholds
    thSet = [repmat({0.5:0.5:40},2,1); {0.4:0.4:20}; spTh];
    % Exclusion thresholds
    sigThSet = [5;5;100;sigTh];
    sMedTh = [1;1;100;sMedTh];
    tMedTh = [3;3;100;tMedTh];
else
    thSet = spTh; sigThSet = sigTh;
end
tmdl = fit_poly([1,Nbt], bvWin, 1);
behTx = ((1:Nbt)'.^[1,0])*tmdl;
% Spontaneous flags
bsFlag = behTx <= 0; brFlag = behTx < brWin;
brFlag = xor(brFlag(:,1),brFlag(:,2));

% Set-point median for the whiskers and nose
% Spontaneous
ssp = [squeeze(median(dlcStack(Nbs:end,bsFlag,:), 2)); zeros(1,Ntr)];
% Trial
tsp = [squeeze(median(dlcStack(Nbs:end,:,:), 2)); zeros(1,Ntr)];
clearvars dlcStack vStack;
% Measuring trials
sSig = cell2mat(cellfun(@(x) squeeze(std(x(bsFlag,:), 1, 1)), behStack, ...
    fnOpts{:}));
sMed = cell2mat(cellfun(@(x) squeeze(median(x(bsFlag,:), 1)), behStack, ...
    fnOpts{:}));
tMed = cell2mat(cellfun(@(x) squeeze(median(x, 1)), behStack, fnOpts{:}));

thrshStrs = arrayfun(@(x,y,z) sprintf("TH s%.2f sp_m%.2f t_m%.2f", x,y,z), ...
    sigThSet, sMedTh, tMedTh);
% Excluding flags: Sigma and median in spontaneous or the whole trial
% median greater than a given threshold
excFlag = sSig > sigThSet | abs(sMed-ssp) > sMedTh | abs(tMed-tsp) > tMedTh;
excFlag = excFlag';

%% Organising figures in subfolders
vwKey = sprintf("VW%.2f - %.2f s", bvWin);
rwKey = sprintf("RW%.2f - %.2f ms", brWin*k);
subFig = "V%.2f - %.2f s R%.2f - %.2f ms";
% Configuration subfolder
subfigDir = fullfile(figureDir, sprintf(subFig, bvWin, brWin*k));
%TODO: Body part subfolder
if exist(subfigDir, "dir")
    figureDir = subfigDir;
else
    if ~mkdir(subfigDir)
        if verbose
            fprintf(1, "Error while creating subfolder!\n")
            fprintf(1, "Placing the figures in 'Figure' directory.\n")
        end
    else
        figureDir = subfigDir;
    end
end
%% Going through the conditions and signals
bfPttrn = "%s %s "+vwKey+" "+rwKey+" EX%d %s";
mbfPttrn = "Mean %s %s"+vwKey+" "+rwKey+" EX%s%s";
mpfPttrn = "Move probability %s %s"+rwKey+" %s";
pfPttrn = "%s %s move probability %.2f "+rwKey+" EX%d %s";
mvdPttrn = "%s dist "+vwKey+" "+rwKey+" EX%s%s";
clMap = lines(Nccond);
% Turning condition flag per trial flags into a page.
pageTrialFlag = reshape(pairedStim, size(pairedStim, 1), 1, []);
% Exclusion flags + paired stimulation/experimental conditions control
% Used trial flags shape #Trial x Behavioural signals x #Conditions
xtf = ~excFlag & pageTrialFlag;
% Excluded flags: #behSign x #Conditions
Nex = reshape(sum(xor(xtf, pageTrialFlag)),Nbs, Nccond);
% Figure names for all signals and conditions
bfNames = arrayfun(@(y) arrayfun(@(x) sprintf(bfPttrn, behNames(x), ...
        consCondNames{y}, Nex(x,y), thrshStrs(x)), 1:Nbs), ...
        1:Nccond, fnOpts{:});
bfNames = cat(1, bfNames{:});
% Figure names for mean trials per condition
mbfNames = arrayfun(@(s) sprintf(mbfPttrn, behNames(s), ...
    sprintf("%s ", consCondNames{:}), sprintf("%d ", Nex(s,:)), ...
    thrshStrs(s)), 1:Nbs);
% Figure names for max value per trial boxchart
mvdNames = arrayfun(@(s) sprintf(mvdPttrn, behNames(s), sprintf("%d ", ...
    Nex(s,:)), thrshStrs(s)), 1:Nbs);
% Computing the SEM and quantiles from the signals per condition.
behSgnls = arrayfun(@(y) arrayfun(@(x) getSEM(behStack{x}, ...
    xtf(:, x, y)), 1:Nbs, fnOpts{:}), 1:Nccond, fnOpts{:});
behSgnls = cat(1, behSgnls{:});
% qSgnls = arrayfun(@(y) arrayfun(@(x) getQs(behStack{x}, xtf(:, x, y)), ...
%     1:Nbs, fnOpts{:}), 1:Nccond, fnOpts{:});
% qSgnls = cat(1, qSgnls{:});
% Getting maximum speed per trial
mvpt = arrayfun(@(y) arrayfun(@(x) ...
    getMaxAbsPerTrial(behStack{x}(:,xtf(:, x, y)), brWin, behTx), ...
    1:Nbs, fnOpts{:}), 1:Nccond, fnOpts{:});
mvpt = cat(1, mvpt{:});
% Crossing thresholds and movement probability
mvFlags = arrayfun(@(y) arrayfun(@(x) compareMaxWithThresh(mvpt{y, x}, ...
    thSet(x)), 1:Nbs, fnOpts{:}), 1:Nccond, fnOpts{:});
mvFlags = cat(1, mvFlags{:}); mov_prob = cellfun(@getAUC, mvFlags);
movSgnls = cellfun(getThreshCross, mvFlags, fnOpts{:});
movSgnls = arrayfun(@(s) cat(1, movSgnls{:,s})', 1:Nbs, fnOpts{:});
ccnMP = arrayfun(@(c) arrayfun(@(s) consCondNames{c}+" "+ ...
    string(mov_prob(c, s)), 1:Nbs), 1:Nccond, fnOpts{:});
ccnMP = cat(1, ccnMP{:});
% Probability figure name with all conditions for all signals
mpfNames = arrayfun(@(s) sprintf(mpfPttrn, behNames(s), ...
    sprintf("%s ", ccnMP(:,s)), thrshStrs(s)), 1:Nbs);
% Trial crossing plot per condition and signals
pfNames = arrayfun(@(y) arrayfun(@(x) sprintf(pfPttrn, behNames(x), ...
        consCondNames{y}, mov_prob(y, x), Nex(x,y), thrshStrs(x)), ...
        1:Nbs), 1:Nccond, fnOpts{:});
pfNames = cat(1, pfNames{:});
resStruct = struct('Name', cellstr(behNames), 'Maximum_value', mvpt, ...
    'MovProb', num2cell(mov_prob));
% Tests for roller movement only 
if Nccond > 1
    cs = 4; prms = nchoosek(1:Nccond,2);
    getDistTravel = @(c) squeeze(sum(abs(behStack{cs}(brFlag, ...
        xtf(:, cs, c))), 1));
    dstTrav = arrayfun(getDistTravel, 1:Nccond, fnOpts{:});
    [pd, hd] = arrayfun(@(p) ranksum(dstTrav{prms(p,1)}, ...
        dstTrav{prms(p,2)}), 1:size(prms,1), fnOpts{:});
    [pm, hm] = arrayfun(@(p) ranksum(mvpt{prms(p,1), cs}, ...
        mvpt{prms(p,2), cs}), 1:size(prms,1), fnOpts{:});
    rsstPttrn = "Results %s %s%s"+rwKey+" EX%s%s.mat";
    rsstName = sprintf(rsstPttrn, behNames(cs), sprintf("%s ", ...
        consCondNames{:}), sprintf("%.2f ", pm{:}), sprintf("%d ", ...
        Nex(cs,:)), thrshStrs(cs));
    rsstPath = behHere(rsstName);
    if ~exist(rsstPath, "file")
        save(rsstPath,"pd", "hd", "pm", "hm");
    end
end

% Plotting results and allocating figure holders
allTrialFigs = gobjects(Nccond, Nbs); muTrialFigs = gobjects(Nbs, 1);
mpFigs = allTrialFigs; mppcFigs = muTrialFigs; bpfFigs = mppcFigs;
for cbs = 1:Nbs
    for ccond = 1:Nccond
        % Each trial and mean overlaid
        figTtl = sprintf("%s %s", behNames(cbs), consCondNames{ccond});
        allTrialFigs(ccond, cbs) = figure('Name', figTtl, 'Color', 'w');
        cax = axes('Parent', allTrialFigs(ccond, cbs), axOpts{:});
        plot(cax, behTx*k, behStack{cbs}(:, xtf(:, cbs, ccond)), ptOpts{1,:}); 
        title(cax, figTtl + " Ex:" + string(Nex(cbs, ccond)));
        mtpObj = plot(cax, behTx*k, behSgnls{ccond, cbs}(:,1), ptOpts{2,:});
        lgObj = legend(mtpObj, behNames(cbs)); set(lgObj, lgOpts{:});
        xlabel(cax, "Time [ms]"); ylabel(cax, yLabels(cbs)); xlim(cax, ...
            bvWin*k)
        mpFigs(ccond, cbs) = plotThetaProgress(mvFlags(ccond, cbs), ...
            thSet(cbs), string(consCondNames{ccond}));
        cax = get(mpFigs(ccond, cbs), "CurrentAxes");
        set(cax, axOpts{:}); xlabel(cax, "\theta "+yLabels(cbs))
        title(cax, sprintf("Trial proportion %s %.3f", figTtl, ...
            mov_prob(ccond, cbs)))
    end
    % Mean for the considered trials
    muTrialFigs(cbs) = figure('Name', "Mean "+behNames(cbs), "Color", "w");
    cax = axes("Parent", muTrialFigs(cbs), axOpts{:});
    arrayfun(@(c) patch(cax, k*behTx([1:end, end:-1:1]),...
        mat2ptch(behSgnls{c, cbs}), 1, phOpts{:}, clMap(c,:)), 1:Nccond)
    pObj = arrayfun(@(c) plot(cax, k*behTx, behSgnls{c, cbs}(:,1), "Color", ...
        clMap(c,:), "LineWidth", 1.5, "DisplayName", consCondNames{c}), ...
        1:Nccond); lgObj = legend(pObj); set(lgObj, lgOpts{:}); xlim(cax, ...
        bvWin*k); xlabel(cax, "Time [ms]"); ylabel(cax, yLabels(cbs)); 
    title(cax, "Mean: "+behNames(cbs))
    % Movement probability for all conditions
    mppcFigs(cbs) = figure('Name', "Move prob "+behNames(cbs), "Color", "w");
    cax = axes("Parent", mppcFigs(cbs), "Colormap", clMap, axOpts{:});
    plot(cax, thSet{cbs}, movSgnls{cbs}, "LineWidth", 0.7)
    xlabel(cax, yLabels(cbs)); ylabel(cax, "Trial crossing"); 
    title(cax, sprintf("Move probability for %s", behNames(cbs)))
    lgObj = legend(ccnMP(:,cbs)); set(lgObj, lgOpts{:});
    ylim(cax, [0, 1])
    % Boxplot for maximum value per condition
    bpfFigs(cbs) = figure('Name', "Maximum value per trial "+behNames(cbs), ...
        'Color',"w"); cax = axes("Parent", bpfFigs(cbs), axOpts{:});
    arrayfun(@(c) boxchart(cax, c*ones(size(mvpt{c,cbs})), ...
        mvpt{c, cbs}, 'Notch', 'on'), 1:Nccond)
    xticks(cax, 1:Nccond); xticklabels(cax, consCondNames)
    mxLevl = ceil(1.05*max(cellfun(@(c) quantile(c, 3/4) + ...
        1.5*iqr(c), mvpt(:, cbs))));
    if mxLevl <= cax.YLim(1)
        mxLevl = cax.YLim(2);
    end
    ylim(cax, [cax.YLim(1), mxLevl]); ylabel(cax, yLabels(cbs));
    title(cax, behNames(cbs)+" max val distribution"); 
end
% Saving the plots
figDir = @(x) fullfile(figureDir, x);
arrayfun(@(c) arrayfun(@(s) saveFigure(allTrialFigs(c, s), ...
    figDir(bfNames(c, s)), 1), 1:Nbs), 1:Nccond);
arrayfun(@(c) arrayfun(@(s) saveFigure(mpFigs(c, s), figDir(pfNames(c, s)), ...
    1), 1:Nbs), 1:Nccond);
arrayfun(@(s) saveFigure(muTrialFigs(s), figDir(mbfNames(s)), 1), 1:Nbs);
arrayfun(@(s) saveFigure(mppcFigs(s), figDir(mpfNames(s)), 1), 1:Nbs);
arrayfun(@(s) saveFigure(bpfFigs(s), figDir(mvdNames(s)), 1), 1:Nbs);

%{
TODO: 
1.- Organise parameters in structures:
    1.1.- Inputs (signal processing, exclusion)
    1.2.- Outputs (Tests, values, maybe stacks)
2.- Organise figures in further subfolders per bodypart
%}
end
%{
function [fig, plObj] = figPlot(x, y, opts, axOpts)
fig = figure; ax = axes('Parent', fig, axOpts{:});
plObj = plot(ax, x, y, opts{:});
end
%}
function [outputArg1,outputArg2] = analyseBehaviour(behDir, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Auxiliary variables
% Input files
dlcPttrn = 'roller*filtered.csv';
rpPttrn = "Roller_position*.csv"; vdPttrn = "roller*.avi";
% Output files
afPttrn = "ArduinoTriggers*.mat"; rfPttrn = "RollerSpeed*.mat";
% DLC variables
dlcVars = {'behDLCSignals', 'dlcNames'};
% lPttrn = "Laser*.csv"; pPttrn = "Puff*.csv";
% Functional cell array
fnOpts = {'UniformOutput', false}; 
% Aesthetic cell arrays
axOpts = {'Box','off','Color','none'};
lgOpts = cat(2, axOpts{1:2}, {'Location','best'});
ptOpts = {"Color", 0.7*ones(1,3), "LineWidth", 0.2;...
    "Color", "k", "LineWidth",  1.5};
% Anonymus functions
behHere = @(x) fullfile(behDir, x);
flfile = @(x) fullfile(x.folder, x.name); iemty = @(x) isempty(x);
istxt = @(x) isstring(x) | ischar(x);
mat2ptch = @(x) [x(1:end,:)*[1;1]; x(end:-1:1,:)*[1;-1]];
getThreshCross = @(x) sum(x)/size(x,1);
getSEM = @(x, idx) [squeeze(mean(x(:,idx),2)), ...
    std(x(:,idx),1,2)./sqrt(sum(idx))];
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
    'MedianTh', tMedTh));

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
    if exist(jointArdPath, "file")
        if verbose
            [~, auxFile] = fileparts(jointArdPath);
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
        atT = arrayfun(@(x, z) cellfun(@(y, a) y+a, ...
            x.atTimes, repmat(z,1,length(x.atTimes)), fnOpts{:}), atV', ...
            num2cell([0, Texp(1:end-1)]), fnOpts{:});
        atT = cat(1, atT{:});
        atTimes = arrayfun(@(x) cat(1, atT{:,x}), 1:size(atT,2), fnOpts{:});
        atNames = atV(1).atNames;
        % Intan
        itT = arrayfun(@(x, z) cellfun(@(y, a) y+a, ...
            x.itTimes, repmat(z,size(x.itTimes)), fnOpts{:}), atV', ...
            num2cell([0, Texp(1:end-1)]), fnOpts{:});
        itT = cat(2, itT{:})';
        itTimes = arrayfun(@(x) cat(1, itT{:,x}), 1:size(itT,2), fnOpts{:});
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
elseif islogical(pairedStim) && size(pairedStim, 1) == Ntr
    % Connected to DE_Jittering
    NTa = sum(pairedStim(:));
    if sum(pairedStim(:)) ~= Ntr
        fprintf(1, "The number of trials in the given flags (%d) do not", ...
            NTa)
        fprintf(1, "correspond to the number of trials in behaviour files");
        fprintf(1, " (%d)!\n", Ntr);
        dbstop in analyseBehaviour at 287
    else
        Nccond = size(pairedStim, 2);
    end
end
% Estimating the setpoint for whiskers and nose signals.
behTx_all = (0:size(behDLCSignals, 1) - 1)/fr;
spFile = fullfile(behDir, "SetPoint.mat");
if ~exist(spFile, "file")
    set_point = arrayfun(@(x) fitSpline(behTx_all, behDLCSignals(:,x), 1, ...
        0.3, 0.95, false), 1:size(behDLCSignals,2), fnOpts{:});
    set_point = [set_point{:}];
    save(spFile, "set_point")
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
behNames = [dlcNames, "Roller speed"];
% Whiskers, nose, and roller speed

%TODO: include user input for the thresholds for whiskers and nose signals.
% Probability thresholds
thSet = [repmat({0.5:0.5:40},2,1); {0.4:0.4:20}; spTh];
% Exclusion thresholds
sigThSet = [5;5;0;sigTh];
sMedTh = [1;1;0;sMedTh];
tMedTh = [1;1;0;tMedTh];

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

% Loop variables
gp = zeros(Nccond, Nbs+1, 'single');
bfPttrn = "%s %s VW%.2f - %.2f s RM%.2f - %.2f ms EX%d %s";
pfPttrn = "%s move probability %.2f RW%.2f - %.2f ms EX%d %s";
behSgnls = cell(Nccond, Nbs); mvFlags = cell(Nccond, Nbs); mvpt = mvFlags;
qSgnls = behSgnls;

% Turning condition flag per trial flags into a page.
pageTrialFlag = reshape(pairedStim, size(pairedStim, 1), 1, []);
% Exclusion flags + paired stimulation/experimental conditions control
% Used trial flags shape #Trial x Behavioural signals x #Conditions
xtf = ~excFlag & pageTrialFlag;
% Excluded flags: #behSign x #Conditions
Nex = squeeze(sum(xor(xtf, pageTrialFlag)));

for ccond = 1:Nccond
    fig = figure('Name',"Trial by trial", "Color","w");
    bfNames = arrayfun(@(x, y, z) sprintf(bfPttrn, consCondNames{ccond}, ...
        x, bvWin, brWin*k, y, z), behNames', Nex(:,ccond), thrshStrs);
    behSgnls{ccond,:} = arrayfun(@(x) getSEM(behStack{x}, xtf(:,x,ccond)), ...
        1:Nbs, fnOpts{:});
end

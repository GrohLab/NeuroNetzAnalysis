function [summStruct, figureDir, behData, aInfo] = ...
    analyseBehaviour(behDir, varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Auxiliary variables
% 'Software' version
softVer = 3.01;
% Input files
dlcPttrn = 'roller*shuffle2*filtered.csv';
% rpPttrn = "Roller_position*.csv"; vdPttrn = "roller*.avi";
% Output files
afPttrn = "ArduinoTriggers*.mat"; rfPttrn = "RollerSpeed*.mat";
% DLC variables
dlcVars = {'behDLCSignals', 'dlcNames', 'vidTx', 'fr', 'refStruct'};
% lPttrn = "Laser*.csv"; pPttrn = "Puff*.csv";
% Functional cell array
fnOpts = {'UniformOutput', false};
% Aesthetic cell arrays
axOpts = {'Box', 'off', 'Color', 'none', 'NextPlot', 'add'};
lgOpts = cat(2, axOpts{1:4}, {'Location', 'best'});
ptOpts = {"Color", 0.7*ones(1,3), "LineWidth", 0.2;...
    "Color", "k", "LineWidth",  1.5};
phOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.25, 'FaceColor'};
% Anonymus functions
behHere = @(x) fullfile(behDir, x);
flfile = @(x) fullfile(x.folder, x.name); iemty = @(x)isempty(x);
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
checkCond = @(x) isscalar(x) && istxt(x) && (strcmpi(x,"P") || strcmpi(x,"L"));

% If the triggers don't represent the same for all the experiment or if you
% want to see only a subset of all the triggers.
defTrigSubs = "all";
checkTrigSubs = @(x) strcmpi(x,"all") || isPositiveIntegerValuedNumeric(x);
% Check paired stimulus
checkPF = @(x) islogical(x) & (ismatrix(x) || isvector(x));
% Check conditions names
checkCondNames = @(x) (iscell(x) || istxt(x)) && numel(x) >= 1;

defSpeedTh = 0.1:0.1:3; % Speed thresholds
checkSpeedTh = @(x) all([diff(x) > 0, x > 0, numel(x) > 1]);

checkLogicFlags = @(x) islogical(x) && isscalar(x);

checkTH = @(x) x > 0 && isscalar(x);
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
addParameter(p, 'verbose', true, checkLogicFlags)
addParameter(p, 'showPlots', true, checkLogicFlags)
addParameter(p, 'figOverWrite', false, checkLogicFlags)

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
spFlag = p.Results.showPlots;
owFlag = p.Results.figOverWrite;

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
    fprintf( 1, 'Time window: %.2f - %.2f ms\n', bvWin*k )
    fprintf( 1, 'Response window: %.2f - %.2f ms\n', brWin*k )
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
    rfName = flfile(rfFiles);

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
    if ~exist('Texp','var')
        if exist('Nt','var')
            Texp = Nt;
        else
            Texp = length(vf)/fr;
        end
    end
end

% DLC signals
%readCSV = @(x) readtable(x, "Delimiter", ",");
dlcFiles = dir( behHere( dlcPttrn ) );
%fid_paths = dir( behHere( "FrameID*.csv" ) );

if ~iemty(dlcFiles)
    behPttrn = "BehaviourSignals%s%s.mat";
    endng = "DLC" + string(extractBetween(dlcFiles(1).name, "DLC", ".csv"));
    behPath = joinBehDates(dlcFiles, "roller", behHere(behPttrn), ...
        endng);

    if ~exist(behPath, "file")
        if verbose
            fprintf(1, "Computing signals... \n")
        end
        % varsInLsr = {'lsrInt', 'delta_tiv', 'Texp_vid', 'Texp_ephys'};
        %{
        % 
    exp_path = getParentDir( behDir, 1);
    tf_paths = dir( fullfile( behDir, "TriggerSignals*.bin") );
    fsf_path = dir( fullfile( behDir, "*_sampling_frequency.mat") );
    if isempty( fsf_path )
        fsf_path = dir( fullfile( exp_path, "ephys*", "*_sampling_frequency.mat") );
    end
    if isempty( tf_paths )
        tf_paths = dir( fullfile( exp_path, "ephys*", "TriggerSignals*.bin") );
    end

    fs_ephys = load( flfile( fsf_path ), "fs" ); fs_ephys = fs_ephys.fs;
    Ns_intan = [tf_paths.bytes]' ./ 4; % 2 signals x 2 bytes per sample.
    Texp_ephys = Ns_intan ./ fs_ephys;
    vidTx = arrayfun(@(x) readCSV( flfile( x ) ), fid_paths, fnOpts{:} );
    vidTx = cellfun(@(x) x.Var2/1e9, vidTx, fnOpts{:} ); % nanoseconds
    if ~isempty(vidTx)
         Texp_vid = cellfun(@(x) diff( x([1,end]) ), vidTx );
           % Difference between intan signals and the video
            delta_tiv = Texp_ephys - Texp_vid;
            delta_tiv( delta_tiv < 0 ) = 0;
        else
            delta_tiv = zeros( numel(dlcFiles), 1 );
        end
        %}

        [lsrInt, delta_tiv, ~, ~, vidTx, trig, dlcTables, fs] = ...
            extractLaserFromVideos( behDir ); %#ok<ASGLU>
        mean_delay = alignVideoWithEphys( lsrInt, trig, fs, behDir );

        unfeas_delay_flag = abs( mean_delay ) > 0.1;
        mean_delay( unfeas_delay_flag ) = delta_tiv( unfeas_delay_flag );
        mean_delay( mean_delay == 0 ) = mean( mean_delay(mean_delay ~= 0) );
        if verbose
            fprintf(1, "Correcting by")
            fprintf(1, " %.2f", mean_delay * k)
            fprintf(1, " ms\n")
        end
        % dlcTables = arrayfun(@(x) readDLCData(flfile(x)), dlcFiles, fnOpts{:});
        [a_bodyParts, refStruct] = cellfun(@(x) getBehaviourSignals(x), ...
            dlcTables, fnOpts{:}); %#ok<ASGLU>
        behStruct = cellfun(@(x) getWhiskANoseFromTable(x), a_bodyParts);
        % Butter order 3 low-pass filter 35 Hz cut off frequency @ 3 dB
        [b, a] = butter(3, 70/fr, 'low');

        % Adding how many frames fell out to the next experiment segment
        behDLCSignals = cell( numel( behStruct ), 1 );
        for x = 1:numel( behStruct )
            if mean_delay(x) > 0
                behDLCSignals{x} = padarray( behStruct(x).Signals, ...
                    round( [mean_delay(x), 0] * fr ), "symmetric", "pre" );
            elseif mean_delay(x) < 0
                behDLCSignals{x} = behStruct(x).Signals( ...
                    round( -mean_delay(x) * fr ):end, : );
            else
                behDLCSignals{x} = behStruct(x).Signals;
            end
            tmDf = length(trig{x})/fs - length(behDLCSignals{x})/fr;
            if tmDf < 0
                behDLCSignals{x} = behDLCSignals{x}(1:end+round( tmDf*fr ),:);
            elseif tmDf > 0
                behDLCSignals{x} = padarray( behDLCSignals{x}, ...
                    round( [tmDf, 0] * fr ), "symmetric", "post" );
            end
        end

        % Concatenating both experiment segments
        % behDLCSignals = cat(1, behStruct(:).Signals);
        behDLCSignals = cat( 1, behDLCSignals{:} );

        % Filter and save the results
        dlcNames = behStruct(1).Names; clearvars behStruct;
        behDLCSignals = filtfilt(b, a, behDLCSignals);

        if ~isreal( behDLCSignals )
            if verbose
                fprintf(1, 'Some signals turned out to be complex. ')
                fprintf(1, 'Transforming to real values...')
            end
            behDLCSignals = real( behDLCSignals );
            if verbose
                fprintf(1, ' Done!\n')
            end
        end

        if verbose
            fprintf(1, 'Saving... ')
        end
        save(behPath, dlcVars{:})
        if verbose
            fprintf(1, "Done!\n")
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
    fr, fr, [], vf*en2cm); [~, Nbt, Ntr] = size(vStack);
% Connection to DE_Jittering or as a stand-alone function.
if istxt(pairedStim) && strcmpi(pairedStim, "none")
    % Stand alone function without stimulus pairing
    pairedStim = true(Ntr, 1);
    Nccond = 1;
    consCondNames = cellstr(atNames(aSub));
elseif islogical(pairedStim)
    % Connected to DE_Jittering
    %NTa = sum(pairedStim(:));
    if size(pairedStim, 1) ~= Ntr
        fprintf(1, "The number of trials in the given flags (%d) do not", ...
            size(pairedStim, 1))
        fprintf(1, " correspond to the number of trials in behaviour files");
        fprintf(1, " (%d)!\n", Ntr);
        fprintf(1, "Removing empty rows... (%d)\n", sum(all(~pairedStim,2)))
        pairedStim(all(~pairedStim,2),:) = [];
    end
    Nccond = size(pairedStim, 2);
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
% Estimating the setpoint for whiskers and nose signals.
behTx_all = (0:size(behDLCSignals, 1) - 1)/fr;
spFile = fullfile(behDir, "SetPoint2.mat");
if ~exist(spFile, "file")
    set_point = cell( Nbs, 1 );
    parfor cs = 1:size( behDLCSignals, 2 )
        set_point{cs} = fitSpline( behTx_all, behDLCSignals( :, cs ), 1, ...
            0.3, 0.95, verbose );
    end
    set_point = [set_point{:}];
    if ~isempty(set_point)
        save(spFile, "set_point")
    end
else
    set_point = load(spFile,"set_point"); set_point = set_point.set_point;
end

% Mean and standard deviation of all signals
[~, mu_dlc, sig_dlc] = zscore( behDLCSignals, 0, 1);
[~, mu_r, sig_r] = zscore( vf*en2cm , 0, 1 );

% Whiskers and nose; video data
[~, dlcStack] = getStacks(false, round(itTimes{iSub}(trigSubs,:)*fr), 'on', ...
    bvWin, fr, fr, [], [behDLCSignals, set_point]');
[~, dlcDiffStack] = getStacks(false, round(itTimes{iSub}(trigSubs,:)*fr), 'on', ...
    bvWin, fr, fr, [], diff(behDLCSignals,1,1)');

behStack = cat(1, dlcStack(1:Nbs,:,:), vStack);
Nbs = size(behStack, 1);
behStack = arrayfun(@(x) squeeze(behStack(x, :, :)), (1:size(behStack,1))',...
    fnOpts{:});
behNames = [dlcNames, "Roller speed"]; behNames(~strlength(behNames))=[];
rollYL = "Roller speed [cm/s]";
if size(behNames,2)>1
    yLabels = [repmat("Angle [°]", 1, Nbs-1), rollYL];
    sym_flag = contains( behNames, "symmetry", "IgnoreCase", true );
    yLabels(sym_flag) = "Symmetry [a.u.]";
else
    yLabels = rollYL;
end

% Movement probability
movStruct = computeMovementProbability(behStack, pairedStim, bvWin, fr);

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
tmdl = fit_poly([1,Nbt], bvWin + [1,-1] * (1/(fr*2)), 1);
behTx = ((1:Nbt)'.^[1,0])*tmdl;

% Spontaneous flags
bsWin = -flip(brWin);
bsFlag = behTx < bsWin; bsFlag = xor(bsFlag(:,1), bsFlag(:,2));

brWin_aux = brWin;
if brWin(1) < 0.12
    brWin_aux = brWin + 0.1;
end

if brWin_aux(2) > bvWin(2)
    if verbose
        fprintf(1, "WARNING! The response window for the stimulated");
        fprintf(1, " whisker gets out of the viewing window!\n");
        fprintf(1, "Reconsider the span of these windows!\n");
    end
    brWin_aux(2) = bvWin(2);
end
brWin = [brWin_aux; repmat( brWin, Nbs-1, 1 )];

brFlag = arrayfun(@(x) behTx < brWin(x,:), 1:Nbs, fnOpts{:} );
brFlag = cellfun(@(x) xor(x(:,1), x(:,2) ), brFlag, fnOpts{:} );
brFlag = cat( 2, brFlag{:} );

% % Set-point median for the whiskers and nose
% % Spontaneous
% ssp = [squeeze(median(dlcStack(Nbs:end,bsFlag,:), 2)); zeros(1,Ntr)];
% % Trial
% tsp = [squeeze(median(dlcStack(Nbs:end,:,:), 2)); zeros(1,Ntr)];
% %clearvars dlcStack vStack;
% % Measuring trials
% sSig = cell2mat(cellfun(@(x) squeeze(std(x(bsFlag,:), 0, 1)), behStack, ...
%     fnOpts{:}));
% sMed = cell2mat(cellfun(@(x) squeeze(median(x(bsFlag,:), 1)), behStack, ...
%     fnOpts{:}));
% tMed = cell2mat(cellfun(@(x) squeeze(median(x, 1)), behStack, fnOpts{:}));

% thrshStrs = arrayfun(@(x,y,z) sprintf("TH s%.2f sp_m%.2f t_m%.2f", x,y,z), ...
%     sigThSet, sMedTh, tMedTh);

% Excluding flags: Sigma and median in spontaneous or the whole trial
% median greater than a given threshold
%excFlag = sSig > sigThSet | abs(sMed-ssp) > sMedTh | abs(tMed-tsp) > tMedTh;
%excFlag = excFlag';
excFlag = false(Ntr, Nbs);

%% Organising figures in subfolders
vwKey = sprintf("V%.2f - %.2f s", bvWin);
rwKey = sprintf("R%.2f - %.2f ms", brWin(2,:)*k);
subFig = "Beh %s %s";
% Configuration subfolder
subfigDir = fullfile(figureDir, sprintf(subFig, vwKey, rwKey));
%TODO: Body part subfolder
metaNameFlag = false;
if exist(subfigDir, "dir")
    figureDir = subfigDir;
else
    if ~mkdir(subfigDir)
        % Print metadata on figure name
        metaNameFlag = true;
        if verbose
            fprintf(1, "Error while creating subfolder!\n")
            fprintf(1, "Placing the figures in 'Figure' directory.\n")
        end
    else
        figureDir = subfigDir;
    end
end
%% Going through the conditions and signals
bfPttrn = "%s %s ";
mbfPttrn = "Mean %s %s";
mpfPttrn = "Move probability %s %s";
pfPttrn = "%s %s move probability %.2f ";
mvdPttrn = "%s dist ";
ttPttrn = "%s in %d from %d trials";
if metaNameFlag
    bfPttrn = bfPttrn + vwKey + " " + rwKey;
    mbfPttrn = mbfPttrn + vwKey + " " + rwKey;
    mpfPttrn = mpfPttrn + rwKey;
    pfPttrn = pfPttrn + rwKey;
    mvdPttrn = mvdPttrn + vwKey + " " + rwKey;
end
bfPttrn = bfPttrn + " EX%d %s";
mbfPttrn = mbfPttrn + " EX%s%s";
mpfPttrn = mpfPttrn + " %s";
pfPttrn = pfPttrn + " EX%d %s";
mvdPttrn = mvdPttrn + " EX%s%s";

% summPttrn = "Summary";
clMap = lines( Nccond );

% Turning condition flag per trial flags into a page.
pageTrialFlag = reshape( pairedStim, size( pairedStim, 1 ), 1, [] );

% Exclusion flags + paired stimulation/experimental conditions control
% Used trial flags shape #Trial x Behavioural signals x #Conditions
xtf = ~excFlag & pageTrialFlag;

% Final number of trials
Na = reshape( sum( xtf, 1 ), Nbs, Nccond );

% Excluded flags: #behSign x #Conditions
Nex = reshape( sum( xor( xtf, pageTrialFlag ) ), Nbs, Nccond );

% Figure names for all signals and conditions
% bfNames = arrayfun( @(y) arrayfun( @(x) sprintf( bfPttrn, behNames(x), ...
%     consCondNames{y}, Nex(x,y), thrshStrs(x)), 1:Nbs ), ...
%     1:Nccond, fnOpts{:} );
bfNames = arrayfun( @(y) arrayfun( @(x) sprintf( bfPttrn, behNames(x), ...
    consCondNames{y}, Nex(x,y) ), 1:Nbs ), 1:Nccond, fnOpts{:} );
bfNames = cat( 1, bfNames{:} );

% Figure names for mean trials per condition
% mbfNames = arrayfun(@(s) sprintf( mbfPttrn, behNames(s), ...
%     sprintf("%s ", consCondNames{:} ), sprintf("%d ", Nex(s,:) ), ...
%     thrshStrs(s) ), 1:Nbs );
mbfNames = arrayfun(@(s) sprintf( mbfPttrn, behNames(s), ...
    sprintf("%s ", consCondNames{:} ), sprintf("%d ", Nex(s,:) ) ), 1:Nbs );

% Figure names for max value per trial boxchart
% mvdNames = arrayfun(@(s) sprintf( mvdPttrn, behNames(s), sprintf("%d ", ...
%     Nex(s,:)), thrshStrs(s) ), 1:Nbs);
mvdNames = arrayfun(@(s) sprintf( mvdPttrn, behNames(s), sprintf("%d ", ...
    Nex(s,:) ) ), 1:Nbs);

% Computing the SEM and quantiles from the signals per condition.
behSgnls = arrayfun(@(y) arrayfun(@(x) getSEM( behStack{x}, ...
    xtf(:, x, y) ), 1:Nbs, fnOpts{:} ), 1:Nccond, fnOpts{:} );
behSgnls = cat( 1, behSgnls{:} );
% qSgnls = arrayfun(@(y) arrayfun(@(x) getQs(behStack{x}, xtf(:, x, y)), ...
%     1:Nbs, fnOpts{:}), 1:Nccond, fnOpts{:});
% qSgnls = cat(1, qSgnls{:});

% Getting maximum speed per trial
mvpt = arrayfun(@(c) arrayfun( @(b) ...
    getMaxAbsPerTrial( behStack{b}(:, xtf(:, b, c)), brWin(b,:), behTx ), ...
    1:Nbs, fnOpts{:} ), 1:Nccond, fnOpts{:} );
mvpt = cat(1, mvpt{:});

% Crossing thresholds and movement probability
ctrl_cond = contains( consCondNames, "Control" ); up_pc = 1.15;
if all( ~ctrl_cond )
    % No control condition found
    mvps = round( max( cellfun( @(x) max( x ), mvpt ) ) * up_pc, 1 );
else
    % Found a control condition
    mvps = round( cellfun( @(x) max( x ), mvpt(ctrl_cond, :) ) * up_pc, 1 );
end
mvps( mvps == 0 ) = 1;

[mvFlags, thSet] = arrayfun( @(y) arrayfun( @(x) ...
    compareMaxWithThresh( mvpt{y, x}, mvps(x) ), 1:Nbs, fnOpts{:} ), ...
    1:Nccond, fnOpts{:} );
thSet = cat( 1, thSet{:} ); thSet = cellfun(@(x) x{:}, thSet(1,:), fnOpts{:} );
mvFlags = cat(1, mvFlags{:}); mov_prob = cellfun( @getAUC, mvFlags );
mov_prob(isnan(mov_prob)) = 0;

% Testing joint amplitude -- CAUTION! Assuming first condition as control!!
[~, mu_ctr, sig_ctr] = cellfun(@(c) zscore( c, 1 ), mvpt(ctrl_cond,:), ...
    fnOpts{:}); mu_ctr = cat(2, mu_ctr{:}); sig_ctr = cat(2, sig_ctr{:});
my_zscore = @(x,m,s) (x - m)./s;
zmvpt = arrayfun(@(b) arrayfun(@(c) ...
    my_zscore( mvpt{c, b}, mu_ctr(b), sig_ctr(b) ), ...
    (1:Nccond)', fnOpts{:} ), 1:Nbs, fnOpts{:} );
zmvpt = cat(2, zmvpt{:});

movSgnls = cellfun(getThreshCross, mvFlags, fnOpts{:});
movSgnls = arrayfun(@(s) cat(1, movSgnls{:,s})', 1:Nbs, fnOpts{:});
ccnMP = arrayfun(@(c) arrayfun(@(s) consCondNames{c}+" "+ ...
    string(mov_prob(c, s)), 1:Nbs), 1:Nccond, fnOpts{:});
ccnMP = cat(1, ccnMP{:});

% Probability figure name with all conditions for all signals
% mpfNames = arrayfun(@(s) sprintf(mpfPttrn, behNames(s), ...
%     sprintf("%s ", ccnMP(:,s)), thrshStrs(s)), 1:Nbs);
mpfNames = arrayfun(@(s) sprintf( mpfPttrn, behNames(s), ...
    sprintf( "%s ", ccnMP(:,s) ) ), 1:Nbs);

% Trial crossing plot per condition and signals
% pfNames = arrayfun(@(y) arrayfun(@(x) sprintf(pfPttrn, behNames(x), ...
%     consCondNames{y}, mov_prob(y, x), Nex(x,y), thrshStrs(x)), ...
%     1:Nbs), 1:Nccond, fnOpts{:});
pfNames = arrayfun(@(y) arrayfun(@(x) sprintf( pfPttrn, behNames(x), ...
    consCondNames{y}, mov_prob(y, x), Nex(x,y) ), ...
    1:Nbs), 1:Nccond, fnOpts{:});
pfNames = cat(1, pfNames{:});

% Auxiliary variables for trial-by-trial plot. Choosing the same number of
% trials per condition for a given behavioural signal
[psSubs, b_c] = find(xtf); psSubs = arrayfun(@(x) psSubs(b_c==x), ...
    unique(b_c), fnOpts{:}); psSubs = reshape(psSubs, Nbs, Nccond);
Nma = min(Na, [], 2); rtSubs = arrayfun(@(b) arrayfun(@(c) ...
    sort(randperm(Na(b,c), Nma(b)),'ascend'), 1:Nccond, fnOpts{:}), ...
    1:Nbs, fnOpts{:}); rtSubs = cat(1, rtSubs{:});
stSubs = cellfun(@(x, y) x(y), psSubs, rtSubs, fnOpts{:});
mvprt = max(cellfun(@(x, y) max(x(y)), mvpt, rtSubs'),[],1);
ttNames = arrayfun(@(b) sprintf(ttPttrn, behNames(b), Nma(b), ...
    Nma(b)+Nex(b)), 1:Nbs);

% Structuring output of the movement
Nae = sum(pairedStim, 1);
movProb = cell2mat(arrayfun(@(ms) sum(ms.MovmentFlags)./Nae, ...
    movStruct, fnOpts{:}));
ms = mat2cell(movStruct, ones(size(movStruct)));
summStruct = cell2mat( arrayfun(@(x) ...
    struct('ConditionName', consCondNames(x), ...
    'Results', struct('BehSigName', cellstr(behNames), ...
    'MaxValuePerTrial', mvpt(x,:), ...
    'Z_Amplitude', zmvpt(x,:), ...
    'AmplitudeIndex', num2cell(mov_prob(x,:)), ...
    'MovProbability', num2cell(movProb(:,x)'),...
    'MovStrucure', ms(:)'), ...
    'NTrials', Nae(x)), 1:Nccond, fnOpts{:}));
summFile = behHere("Simple summary.mat");

% Exporting data
behData = struct( 'Data', cat( 3, behStack{:} ), 'Conditions', pairedStim, ...
    'ConditionNames', string(consCondNames), 'Considered', xtf, ...
    'ZWhole', [mu_dlc(:), sig_dlc(:); mu_r(:), sig_r(:)], ...
    'Zms', [mu_ctr; sig_ctr] );

aInfo = struct('AnalysisTime', datetime(), 'Version', softVer, ...
    'Evoked', rwKey, 'VieWin', vwKey);
if exist(summFile,"file")
    behFile = load(summFile);
    if ~any(contains(fieldnames(behFile), 'aInfo')) || ...
            ~compareSoftVers(behFile.aInfo, aInfo)
        save(summFile, "summStruct", "aInfo", "-append")
    else
        fprintf(1, "File exists, not saving summary numbers.\n")
    end
else
    save(summFile, "summStruct", "aInfo", "-v7.3")
end
% Tests for roller movement only
if Nccond > 1
    cs = contains( behNames, "Roller" );
    prms = nchoosek(1:Nccond,2);
    getDistTravel = @(c) squeeze(sum(abs(behStack{cs}(brFlag(:, cs), ...
        xtf(:, cs, c))), 1));
    dstTrav = arrayfun( getDistTravel, 1:Nccond, fnOpts{:} );
    if all(~cellfun(@isempty, dstTrav),"all")
        [pd, hd] = arrayfun(@(p) ranksum(dstTrav{prms(p,1)}, ...
            dstTrav{prms(p,2)}), 1:size(prms,1), fnOpts{:});
        [pm, hm] = arrayfun(@(p) ranksum(mvpt{prms(p,1), cs}, ...
            mvpt{prms(p,2), cs}), 1:size(prms,1), fnOpts{:});
        rsstPttrn = "Results %s %s%s" + rwKey + " EX%s%s.mat";
        % rsstName = sprintf( rsstPttrn, behNames(cs), sprintf("%s ", ...
        %     consCondNames{:} ), sprintf( "%.2f ", pm{:} ), sprintf("%d ", ...
        %     Nex(cs,:) ), thrshStrs(cs) );
        rsstName = sprintf( rsstPttrn, behNames(cs), sprintf("%s ", ...
            consCondNames{:} ), sprintf( "%.2f ", pm{:} ), sprintf("%d ", ...
            Nex(cs,:) ) );
        rsstPath = behHere(rsstName);
        if ~exist(rsstPath, "file")
            try
                save( rsstPath,"pd", "hd", "pm", "hm" );
            catch
                % rsstName = sprintf( rsstPttrn, behNames(cs), sprintf("%d ", ...
                %     numel(consCondNames) ), sprintf( "%.2f ", pm{:}), sprintf("%d ", ...
                %     Nex(cs,:)), thrshStrs(cs) );
                rsstName = sprintf( rsstPttrn, behNames(cs), ...
                    sprintf( "%d ", numel(consCondNames) ), ...
                    sprintf( "%.2f ", pm{:}), sprintf( "%d ", Nex(cs,:) ) );
                rsstPath = behHere( rsstName );
                save( rsstPath,"pd", "hd", "pm", "hm" );
            end
        end
    end
end

%% Figure creation
% Plotting results and allocating figure holders
allTrialFigs = gobjects(Nccond, Nbs); muTrialFigs = gobjects(Nbs, 1);
mpFigs = allTrialFigs; mppcFigs = muTrialFigs; bpfFigs = mppcFigs;
tptFigs = gobjects(Nbs,1); ttax = gobjects(Nccond,1);
createtiles = @(f, r, c) tiledlayout( f, r, c, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');
for cbs = 1:Nbs
    tptFigs(cbs) = figure('Name', behNames(cbs), 'Color', 'w', ...
        'Visible', spFlag);
    ttpt = createtiles( tptFigs(cbs), 1, Nccond );
    for ccond = 1:Nccond
        % Equal (random) trials for all conditions
        ttax(ccond) = nexttile(ttpt); set( ttax(ccond), axOpts{:} )
        %subplot(1, Nccond, ccond, axOpts{:}, ...
        %'Parent', tptFigs(cbs));
        imagesc(ttax(ccond), behTx*k, [], ...
            behStack{cbs}(:,stSubs{cbs, ccond})'/mvprt(cbs) )
        title(ttax(ccond), consCondNames{ccond});
        ylim(ttax(ccond), [0.5,Nma(cbs)+0.5]);
        xlim(ttax(ccond), behTx([1,end])*k);
        colormap( ttax(ccond), bwg_cm(256) )

        % Each trial and mean overlaid
        figTtl = sprintf("%s %s", behNames(cbs), consCondNames{ccond});
        allTrialFigs(ccond, cbs) = figure('Name', figTtl, 'Color', 'w',...
            'Visible', spFlag);
        tatf = createtiles( allTrialFigs(ccond, cbs), 1, 1);
        % cax = axes('Parent', allTrialFigs(ccond, cbs), axOpts{:});
        cax = nexttile( tatf ); set( cax, axOpts{:} )
        plot(cax, behTx*k, behStack{cbs}(:, xtf(:, cbs, ccond)), ptOpts{1,:});
        title(cax, figTtl + " Ex:" + string(Nex(cbs, ccond)));
        mtpObj = plot(cax, behTx*k, behSgnls{ccond, cbs}(:,1), ptOpts{2,:});
        lgObj = legend(mtpObj, behNames(cbs)); set(lgObj, lgOpts{:});
        xlabel(cax, "Time [ms]"); ylabel(cax, yLabels(cbs)); xlim(cax, ...
            bvWin*k)
        mpFigs(ccond, cbs) = plotThetaProgress(mvFlags(ccond, cbs), ...
            thSet(cbs), string(consCondNames{ccond}), 'showPlots', spFlag);
        cax = get(mpFigs(ccond, cbs), "CurrentAxes");
        set(cax, axOpts{:}); xlabel(cax, "\theta "+yLabels(cbs))
        title(cax, sprintf("Trial proportion %s %.3f", figTtl, ...
            mov_prob(ccond, cbs)))
    end
    xlabel(ttax, "Time [ms]"); ylabel(ttax(1), "Trials");
    ttax(Nccond).YAxis.Visible = 'off';
    colorbar(ttax(Nccond), axOpts{1:2}, 'Location', 'west', ...
        'TickDirection', 'none', 'Ticks', [13, 243], 'TickLabels', ...
        {'Backward', 'Forward'}, 'AxisLocation', 'out');

    % Mean for the considered trials
    muTrialFigs(cbs) = figure('Name', "Mean "+behNames(cbs), ...
        "Color", "w", 'Visible', spFlag);
    tmtf = createtiles( muTrialFigs(cbs), 1, 1);
    %cax = axes("Parent", muTrialFigs(cbs), axOpts{:});
    cax = nexttile( tmtf ); set( cax, axOpts{:} );
    arrayfun(@(c) patch(cax, k*behTx([1:end, end:-1:1]),...
        mat2ptch(behSgnls{c, cbs}), 1, phOpts{:}, clMap(c,:)), 1:Nccond)
    pObj = arrayfun(@(c) plot(cax, k*behTx, behSgnls{c, cbs}(:,1), "Color", ...
        clMap(c,:), "LineWidth", 1.5, "DisplayName", consCondNames{c}), ...
        1:Nccond); lgObj = legend(pObj); set(lgObj, lgOpts{:}); xlim(cax, ...
        bvWin*k); xlabel(cax, "Time [ms]"); ylabel(cax, yLabels(cbs));
    title(cax, "Mean: "+behNames(cbs))

    % Movement probability for all conditions
    mppcFigs(cbs) = figure('Name', "Move prob "+behNames(cbs), ...
        "Color", "w", 'Visible', spFlag);
    tmppc = createtiles( mppcFigs(cbs), 1, 1 );
    % cax = axes("Parent", mppcFigs(cbs), "Colormap", clMap, axOpts{:});
    cax = nexttile(tmppc); set( cax, "Colormap", clMap, axOpts{:});
    plot(cax, thSet{cbs}, movSgnls{cbs}, "LineWidth", 0.7);
    xlabel(cax, yLabels(cbs)); ylabel(cax, "Trial crossing");
    title(cax, sprintf("Move probability for %s", behNames(cbs)))
    lgObj = legend(ccnMP(:,cbs)); set(lgObj, lgOpts{:});
    axis( cax, [thSet{cbs}([1, end]); 0; 1] )

    % Boxplot for maximum value per condition
    bpfFigs(cbs) = figure('Name', ...
        join( ["Maximum value per trial", behNames(cbs)] ), ...
        'Color', "w", 'Visible', spFlag);
    tbpf = createtiles( bpfFigs(cbs), 1, 1 );
    %cax = axes("Parent", bpfFigs(cbs), axOpts{:} );
    cax = nexttile(tbpf); set( cax, axOpts{:} );
    arrayfun(@(c) boxchart(cax, c*ones(size(mvpt{c,cbs})), ...
        mvpt{c, cbs}, 'Notch', 'on'), 1:Nccond)
    xticks(cax, 1:Nccond); xticklabels(cax, consCondNames)
    mxLevl = ceil(1.05*max(cellfun(@(c) quantile(c, 3/4) + ...
        1.5*iqr(c), mvpt(:, cbs))));
    if mxLevl <= cax.YLim(1)
        mxLevl = cax.YLim(2);
    end
    ylim(cax, [cax.YLim(1), mxLevl]); ylabel(cax, yLabels(cbs));
    title(cax, join( [behNames(cbs), "L\infty", ...
        sprintf("%.2f - %.2f ms", brWin(cbs,:) * k )] ) );
end
% Saving the plots
figDir = @(x) fullfile(figureDir, x);

arrayfun(@(c) arrayfun(@(s) saveFigure(allTrialFigs(c, s), ...
    figDir(bfNames(c, s)), 1, owFlag ), 1:Nbs), 1:Nccond);
arrayfun(@(c) arrayfun(@(s) saveFigure(mpFigs(c, s), figDir(pfNames(c, s)), ...
    1, owFlag), 1:Nbs), 1:Nccond);
arrayfun(@(s) saveFigure(muTrialFigs(s), figDir(mbfNames(s)), 1, owFlag ), 1:Nbs);
arrayfun(@(s) saveFigure(mppcFigs(s), figDir(mpfNames(s)), 1, owFlag ), 1:Nbs);
arrayfun(@(s) saveFigure(bpfFigs(s), figDir(mvdNames(s)), 1, owFlag ), 1:Nbs);
arrayfun(@(s) saveFigure(tptFigs(s), figDir(ttNames(s)), 1, owFlag ), 1:Nbs);
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

i = 1;
for t = find(pairedStim(:,1))'
figure(i); plot(behTx, behStack{4}(:,t), "Color", 0.4*ones(1,3)); hold on; plot(behTx(behTx<0), behTx(behTx<0).^[1,0] * smdl(:,i), 'b');
plot(behTx(behTx>0), behTx(behTx>0).^[1,0] * emdl(:,i), "r")
ssme = sqrt(sum(((behStack{4}(behTx<0,t) - (behTx(behTx<0).^[1,0] * smdl(:,i))).^2)))/abs(behTx(1));
esme = sqrt(sum(((behStack{4}(behTx>0,t) - (behTx(behTx>0).^[1,0] * emdl(:,i))).^2)))/abs(behTx(end));
xline(0, '--k')
[~, mxVel] = min(abs(zipt{t,2}));
xline(zipt{t,2}(mxVel),'g')
scatter(zipt{t,2}, interp1(behTx, behStack{4}(:,t), zipt{t,2}), '.m')
scatter(zipt{t,1}, interp1(behTx, behStack{4}(:,t), zipt{t,1}), 'xk')
text(behTx(1), behStack{4}(1,t), string(ssme))
text(behTx(end), behStack{4}(end,t), string(esme))
title(sprintf('%d, %.3f, %.3f', rsm(t,1), mvpt_rs(t), (-ssme+esme)/(esme+ssme)))
i = i+1;
end

%}